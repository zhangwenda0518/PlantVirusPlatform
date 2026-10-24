# -*- coding: utf-8 -*-
"""迁移弧线动画 → GIF（plotly + kaleido + Pillow，全离线）。

用途（用户点名要 phymapr 的「GIF 动画经纬度地图」）：网页里的时间滑杆动画
（tools.html `drawPgAnim`）只能在浏览器里看，导出成 GIF 才能进 PPT/补充材料。
本模块在**服务端**把同一批数据渲染成 GIF：

  · 帧 = MOTP 的时间窗（summary.motp 的 bins/rows，与 motp.tsv 同口径）；
  · 每帧画该时间窗内的迁移弧线 —— 弧线几何与网页 `_pgArc` **同一套**
    （两端点 + 中点法向偏移控制点的二次贝塞尔）；
  · 弧线颜色 = 迁移事件分类的走廊主导档（与网页地图同色）；
    线宽按**全程**最大计数归一 —— 逐帧归一会让安静窗口里的少量迁移
    看起来和暴发期一样粗，动画反而骗人；
  · 区域锚点（点大小 ∝ 序列数）与标签每帧都在，读者能定位起点终点。

离线：plotly 的 world 底图是内置 TopoJSON、kaleido 本地渲染、Pillow 拼帧，
全程零网络请求（与平台的离线约束一致）。

依赖 plotly / kaleido / Pillow —— 都是环境里已有的（不新增安装项）；
在函数内部 import，常规分析路径不为它付启动成本。

不静默：帧数上限、被跳过的无坐标弧线、无 MOTP 数据的运行，全部进返回
dict（n_frames / frames_dropped / n_arcs_skipped / note），前端如实展示。
"""
import io
import math
import os
import threading
import time

# 弧线几何与前端 _pgArc 保持一致（sag 越大弧越鼓）
_ARC_SAG = 0.18
_ARC_SEG = 24
# 迁移事件分类的走廊主导档配色（与 tools.html CLS_C 同色）
_CLASS_RGB = {'Direct': '26,127,55', 'Indirect': '185,119,14',
              'Distant': '125,60,152', 'Unresolved': '127,140,141'}
_CLASS_CN = {'Direct': 'Direct', 'Indirect': 'Indirect', 'Distant': 'Distant',
             'Unresolved': 'No coords'}
_DEFAULT_RGB = '192,57,43'          # 老 run 没有分类数据时的单色红
_GEO_CONF = dict(scope='world', projection_type='natural earth',
                 showcountries=True, countrycolor='#d5dbdb',
                 showland=True, landcolor='#f6f7f9',
                 showocean=True, oceancolor='#eaf2fb',
                 showcoastlines=True, coastlinecolor='#b8c4cc')


def arc_curve(a, b, seg=_ARC_SEG, sag=_ARC_SAG):
    """两端点 [lat, lon] → 二次贝塞尔弧线的 (lats, lons)。

    与前端 `_pgArc` 同一几何：控制点在两端中点、沿连线法向偏移 sag 比例。
    网页/GIF/后续导出共用一套弧线，避免同一个迁移在两张图里走不同路线。

    ⚠️ 2026-10-03 修：经度插值走**最短方向** —— 两端经度差 >180° 时把西端
    +360（差 <−180° 时把东端 −360）再插值。原来直接用原始经度插值，
    178 → −179.5 这种跨 ±180° 日界线的锚点会沿"绕地球一圈"的方向穿过经度 0，
    弧线横穿整张地图（而两端其实只隔 2.5°）。插值后落在 [−180, 180] 外、
    等价于绕回另一侧的点由 plotly 按等价经度绘制，落点不变。
    """
    alat, alon = float(a[0]), float(a[1])
    blat, blon = float(b[0]), float(b[1])
    dlon = blon - alon
    if dlon > 180.0:            # 最短方向：把 b 挪到 a 的西侧
        blon -= 360.0
    elif dlon < -180.0:         # 最短方向：把 b 挪到 a 的东侧
        blon += 360.0
    mlat, mlon = (alat + blat) / 2.0, (alon + blon) / 2.0
    clat = mlat - sag * (blon - alon)
    clon = mlon + sag * (blat - alat)
    lats, lons = [], []
    for i in range(seg + 1):
        t = i / float(seg)
        u = 1.0 - t
        lats.append(u * u * alat + 2 * u * t * clat + t * t * blat)
        lons.append(u * u * alon + 2 * u * t * clon + t * t * blon)
    return lats, lons


def _bin_key(v):
    """分箱左界 → 对齐键（浮点分箱在 JSON 里会带尾差，先取整到毫秒级）。"""
    return round(float(v) * 1000.0)


def _axis_text(axis):
    return {'dated_tree': '定年树节点年代',
            'tip_years': '后代叶采样年（近似）'}.get(axis, '时间轴未知')


def _clamp_int(v, lo, hi, default):
    try:
        v = int(v)
    except (TypeError, ValueError):
        return default
    return max(lo, min(v, hi))


# kaleido 的浏览器是进程级全局单例：不启动 server 时**每次** to_image 都要
# 冷启一个浏览器（实测 4~5 秒/帧，与分辨率无关），启动后首帧 4 秒、
# 其后每帧约 0.2 秒。server 只能开一次/关一次，故用引用计数：并发导出
# （工具页允许）共享它，最后一个退出者负责关；异常路径靠 with 兜底。
_KA_LOCK = threading.Lock()
_KA_USERS = 0


class _kaleido_server(object):
    """引用计数的 kaleido 常驻渲染进程（拿不到就算了，退回逐帧冷启）。"""

    warm = False

    def __enter__(self):
        global _KA_USERS
        try:
            import kaleido
        except ImportError:
            return self
        with _KA_LOCK:
            if _KA_USERS == 0:
                try:
                    kaleido.start_sync_server()
                except Exception:       # noqa: BLE001 —— 降级：仍能逐帧渲染
                    return self
            _KA_USERS += 1
        self.warm = True
        return self

    def __exit__(self, *exc):
        global _KA_USERS
        if not self.warm:
            return False
        with _KA_LOCK:
            _KA_USERS = max(0, _KA_USERS - 1)
            if _KA_USERS == 0:
                try:
                    import kaleido
                    kaleido.stop_sync_server()
                except Exception:       # noqa: BLE001
                    pass
        return False


def render_available():
    """(ok, reason)：本进程能否离线渲染 GIF（plotly + kaleido + Pillow）。

    打包版（dev_tools/VirusPlatform.spec）默认不含 kaleido——按钮点下去要
    给一句人话，而不是 ImportError 堆栈。kaleido 找 Chrome 的失败只能在
    真正渲染时暴露，这里管不到，由 build_migration_gif 首帧兜底提示。
    """
    try:
        import plotly.graph_objects      # noqa: F401
    except ImportError as e:
        return False, f'服务端缺少 plotly，无法导出 GIF（{e}）'
    try:
        import PIL                       # noqa: F401
    except ImportError as e:
        return False, f'服务端缺少 Pillow，无法拼帧导出 GIF（{e}）'
    try:
        import kaleido                   # noqa: F401
    except ImportError as e:
        return False, ('服务端缺少 kaleido，无法渲染 GIF；'
                       f'源码版执行 pip install kaleido 即可（{e}）')
    return True, ''


def plan_frames(summary, max_frames=80):
    """summary → 逐帧数据（不渲染）：[{lo, hi, tot, arcs[]}] + 统计。

    抽出来单独可用：测试与前端「预计帧数」都靠它，渲染只是最后一步。
    弧线按计数降序，超过 max_arcs 的截断并计数（不静默）。
    """
    m = summary.get('motp') or {}
    bins = m.get('bins') or []
    rows = m.get('rows') or []
    coords = summary.get('coords') or {}
    if not bins or not rows:
        raise ValueError(
            '这次运行没有 MOTP 分箱数据（motp_bin=0 或没有可计年的迁移）——'
            '先勾上「迁移随时间（MOTP）」并设分箱宽度后重跑，再导出 GIF。')
    usable = [k for k in coords if k]
    if len(usable) < 2:
        miss = summary.get('coord_missing') or []
        raise ValueError(
            f'可用坐标区划只有 {len(usable)} 个，连不成弧线'
            + (f'（缺坐标：{"、".join(miss)}）' if miss else '')
            + '；请在上方「区域坐标表」补坐标后重跑。')
    by_bin = {_bin_key(b['lo']): [] for b in bins}
    n_skip = 0
    n_offbin = 0
    for r in rows:
        k = _bin_key(r.get('lo'))
        if k not in by_bin:
            # 正常不该发生（rows 与 bins 同源）；真发生了要计数上报，
            # 否则「迁移总数对不上」会被当成数据本身的问题
            n_offbin += 1
            continue
        if r.get('from') not in coords or r.get('to') not in coords:
            n_skip += 1
            continue
        by_bin[k].append(r)
    n_bins_all = len(bins)
    step = max(1, int(math.ceil(n_bins_all / float(max_frames))))
    picked = bins[::step]
    frames = []
    for b in picked:
        rs = by_bin[_bin_key(b['lo'])]
        frames.append({'lo': b['lo'], 'hi': b['hi'],
                       'tot': sum(float(r.get('count') or 0) for r in rs),
                       'arcs': sorted(rs,
                                      key=lambda r: -(float(r.get('count')
                                                            or 0)))})
    return {
        'frames': frames,
        'n_bins': n_bins_all,
        'frames_dropped': n_bins_all - len(picked),
        'step': step,
        'n_arcs_skipped': n_skip,
        'n_rows_offbin': n_offbin,
        'axis': m.get('axis'),
        'bin_width': m.get('bin_width'),
        'note': m.get('note') or '',
        'n_timed': m.get('n_timed'),
        'n_transitions': m.get('n_transitions'),
    }


def _figure(summary, fr, plan, width, height, max_arcs):
    """单帧的 plotly Figure（区域锚点 + 本窗弧线 + 文案）。"""
    import plotly.graph_objects as go

    coords = summary.get('coords') or {}
    n_of = summary.get('regions') or {}
    cls_of = {(c.get('from'), c.get('to')): c
              for c in (summary.get('corridor_classes') or [])}
    # 线宽按全程最大计数归一（帧间可比：安静窗口的少量迁移就是细线）
    vmax = 1.0
    for f2 in plan['frames']:
        for r in f2['arcs']:
            vmax = max(vmax, float(r.get('count') or 0))

    fig = go.Figure()
    names = sorted(coords)
    mx = max([float(n_of.get(k) or 0) for k in names] + [1.0])
    fig.add_trace(go.Scattergeo(
        lon=[coords[k][1] for k in names], lat=[coords[k][0] for k in names],
        mode='markers+text', text=names, textposition='top center',
        textfont=dict(size=11),
        marker=dict(size=[7 + math.sqrt(float(n_of.get(k) or 0) / mx) * 17
                          for k in names],
                    color='#2e86c1', line=dict(color='#ffffff', width=1)),
        hoverinfo='skip', showlegend=False))

    arcs = fr['arcs'][:max_arcs]
    tx, ty, tt = [], [], []
    # scattergeo 的 line.color/width 不接受逐点数组 —— 一条弧一个 trace
    # （每帧 ≤ max_arcs 条，开销可忽略）；按计数升序画，粗弧压在细弧上面
    for r in sorted(arcs, key=lambda x: float(x.get('count') or 0)):
        a, b = coords[r['from']], coords[r['to']]
        la, lo = arc_curve(a, b)
        cl = (cls_of.get((r['from'], r['to'])) or {}).get('label')
        fig.add_trace(go.Scattergeo(
            lon=lo, lat=la, mode='lines', hoverinfo='skip', showlegend=False,
            line=dict(color='rgba(%s,0.85)' % _CLASS_RGB.get(cl, _DEFAULT_RGB),
                      width=1.0 + 7.0 * (float(r.get('count') or 0) / vmax))))
    # 标签：弧线少的时候画两端名+次数，多了会糊成一团（>6 条只画细弧）。
    # 落点在弧顶再沿隆起方向外推 0.6 倍矢高 —— 贝塞尔弧顶就是离弦最远点，
    # 直接放弧顶会压在弧线上（多帧动画里粗弧线宽可达 8px，标签必然被盖）。
    if arcs and len(arcs) <= 6:
        for r in arcs:
            a, b = coords[r['from']], coords[r['to']]
            la, lo = arc_curve(a, b)
            mi = len(la) // 2
            tx.append(lo[mi] + (lo[mi] - (lo[0] + lo[-1]) / 2.0) * 0.6)
            ty.append(la[mi] + (la[mi] - (la[0] + la[-1]) / 2.0) * 0.6)
            tt.append('%s→%s %s' % (r['from'], r['to'], r.get('count')))
        fig.add_trace(go.Scattergeo(
            lon=tx, lat=ty, mode='text', text=tt,
            textfont=dict(size=10, color='#333333'), hoverinfo='skip',
            showlegend=False))

    # 弧线颜色图例（只列本图实际出现的档）；plotly annotation 支持
    # <span style="color:...">，与网页地图的图例文案保持一致
    seen = []
    for r in fr['arcs']:
        cl = (cls_of.get((r['from'], r['to'])) or {}).get('label')
        if cl and cl not in seen:
            seen.append(cl)
    legend = '　'.join('<span style="color:#%s">■</span> %s'
                       % ({'Direct': '1a7f37', 'Indirect': 'b9770e',
                           'Distant': '7d3c98', 'Unresolved': '7f8c8d'}
                          .get(k, 'c0392b'), _CLASS_CN.get(k, k))
                       for k in seen) if seen else ''

    ann = [dict(x=0.01, y=1.06, xref='paper', yref='paper', showarrow=False,
                align='left', font=dict(size=12),
                text='Window <b>%s - %s</b>　migrations <b>%g</b>%s'
                     % (fr['lo'], fr['hi'], fr['tot'],
                        ('(%d more not drawn)' % (len(fr['arcs']) - max_arcs)
                         if len(fr['arcs']) > max_arcs else '')))]
    if legend:
        ann.append(dict(x=1.0, y=1.06, xref='paper', yref='paper',
                        showarrow=False, align='right', font=dict(size=11),
                        text=legend))
    n_dup = summary.get('n_sample_coords') or 0
    foot = ('Time axis: %s　bin %s y　%d windows (1 frame per %d)'
            % (_axis_text(plan['axis']), plan['bin_width'], plan['n_bins'],
               plan['step'])
            + ('　%d sample coords' % n_dup if n_dup else ''))
    if plan['n_arcs_skipped']:
        foot += '　<b>%d migrations missing coords, not drawn</b>' % plan['n_arcs_skipped']
    if plan.get('n_rows_offbin'):
        foot += ('　<b>%d migrations off-bin, not drawn</b>'
                 % plan['n_rows_offbin'])
    ann.append(dict(x=0.01, y=-0.04, xref='paper', yref='paper',
                    showarrow=False, align='left', font=dict(size=10,
                                                             color='#666666'),
                    text=foot))

    fig.update_layout(width=width, height=height, showlegend=False,
                      margin=dict(l=0, r=0, t=34, b=26),
                      font=dict(size=11), geo=dict(_GEO_CONF), annotations=ann)
    # 自动视野：点/弧坐标跨度远小于全球（单一国家/省域内迁移）时，
    # 把底图裁到数据范围并留 ~18% 边距；跨洲数据仍用世界视野
    lons, lats = [], []
    for tr in fig.data:
        lo = getattr(tr, "lon", None)
        la = getattr(tr, "lat", None)
        if lo is not None and la is not None:
            lons += [v for v in lo if v is not None]
            lats += [v for v in la if v is not None]
    if lons and lats and (max(lons) - min(lons)) < 120             and (max(lats) - min(lats)) < 90:
        pad_lo = max(3.0, (max(lons) - min(lons)) * 0.18)
        pad_la = max(3.0, (max(lats) - min(lats)) * 0.18)
        fig.layout.geo.lonaxis.range = [min(lons) - pad_lo, max(lons) + pad_lo]
        fig.layout.geo.lataxis.range = [min(lats) - pad_la, max(lats) + pad_la]
    return fig


def build_migration_gif(summary, out_path, width=960, height=520, fps=2,
                        max_arcs=24, max_frames=80, progress=None,
                        cancel=None):
    """summary(dict) + 输出路径 → 写 GIF；返回统计 dict（含原因说明）。

    fps      帧率：每帧停留 1000/fps 毫秒（默认 2 → 每帧 0.5 秒）
    max_arcs 单帧最多画多少条弧线（多了线会糊；超出条数如实写进帧标题）
    max_frames GIF 帧数上限：窗数多于它时按等间隔取样，丢了几窗如实上报
    progress 可选回调 (frac, msg)（0–1；渲染是逐帧的，用它喂日志）
    cancel   可选 threading.Event：帧间检查，置位即抛错终止（临时文件不落地）
    """
    from PIL import Image

    ok, why = render_available()
    if not ok:
        raise ValueError(why)
    plan = plan_frames(summary, max_frames=max_frames)
    frames = []
    with _kaleido_server() as ka:
        for i, fr in enumerate(plan['frames']):
            if cancel is not None and cancel.is_set():
                raise RuntimeError('导出已取消')
            fig = _figure(summary, fr, plan, width, height, max_arcs)
            try:
                png = fig.to_image(format='png')
            except Exception as e:      # noqa: BLE001
                if i == 0:
                    raise RuntimeError(
                        'GIF 渲染失败：%s（若为找不到浏览器，请安装 Chrome '
                        '或 Edge，或改用源码版环境）' % e) from e
                raise
            frames.append(Image.open(io.BytesIO(png)).convert('RGB'))
            if progress:
                progress((i + 1) / float(len(plan['frames'])),
                         '渲染 GIF 第 %d/%d 帧（时间窗 %s–%s，本窗迁移 %g 次）'
                         % (i + 1, len(plan['frames']), fr['lo'], fr['hi'],
                            fr['tot']))
    if not frames:
        raise ValueError('没有可渲染的帧（MOTP 分箱为空）')
    # 先写 .part 再原子替换：旧 GIF 在导出完成前始终可下载，中途失败/取消
    # 也不会留下半截文件被当成品供应。
    # ⚠️ 2026-10-03：整段包 try/finally —— 渲染/保存/替换任一步抛错（含
    # 取消）都**不留 .part 残留**。原来失败时 .part 会一直躺在输出目录里，
    # 下次导出看不到"上次是怎么坏的"，还占磁盘。
    tmp_path = out_path + '.part'
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or '.',
                exist_ok=True)
    try:
        frames[0].save(tmp_path, format='GIF', save_all=True,
                       append_images=frames[1:],
                       duration=max(40, int(round(1000.0 / max(0.5, fps)))),
                       loop=0, optimize=True)
        # Windows 上 os.replace 会被「另一个进程正打开着目标文件」挡下
        # （PermissionError）—— 典型场景：上一次导出的 GIF 正在被浏览器下载。
        # 短暂重试足够跨过下载尾段；持续被占则把错误如实抛出（别吞）。
        for _try in range(5):
            try:
                os.replace(tmp_path, out_path)
                break
            except PermissionError:
                if _try == 4:
                    raise
                time.sleep(0.4)
    finally:
        # 成功时 os.replace 已把 .part 移走（此处不存在）；失败/取消时清掉残留
        if os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
    return {
        'path': out_path,
        'file': os.path.basename(out_path),
        'bytes': os.path.getsize(out_path),
        'n_frames': len(frames),
        'n_bins': plan['n_bins'],
        'frames_dropped': plan['frames_dropped'],
        'every_n_bins': plan['step'],
        'n_arcs_skipped': plan['n_arcs_skipped'],
        'n_rows_offbin': plan['n_rows_offbin'],
        'axis': plan['axis'],
        'bin_width': plan['bin_width'],
        'n_timed': plan['n_timed'],
        'n_transitions': plan['n_transitions'],
        'width': width, 'height': height, 'fps': fps,
        'warm_render': ka.warm,
        'note': plan['note'],
        'method': ('帧 = MOTP 时间窗；弧线 = 该窗内的迁移计数（颜色 = 迁移事件'
                   '分类的走廊主导档，线宽按全程最大计数归一）；'
                   '底图为 plotly 内置 TopoJSON，kaleido 本地渲染（离线）。'),
    }
