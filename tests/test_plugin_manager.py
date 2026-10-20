# -*- coding: utf-8 -*-
"""插件系统单测：安装 / 升级 / 卸载 / zip 路径穿越防护 / Cfg 填缺合并。

plugins/ 与安装根全部落在 pytest tmp_path，不触碰真实 plugins/ 目录；
Cfg 合并逻辑用 Config.__new__ 跳过重初始化，只挂被测属性。
"""
import json
import os
import zipfile

import pytest

from Virus_Platform_Core import plugin_manager as pm


# ------------------------------------------------------------------
# 造包工具
# ------------------------------------------------------------------
def _make_plugin_zip(path, name='demo-tool', version='1.0', kind='tool',
                     mapping=None, layout='subdir', manifest_extra=None,
                     extra_member=None):
    """构造插件 zip。layout: subdir=唯一顶层目录 / flat=清单在顶层。"""
    top = name if layout == 'subdir' else ''
    body = {
        'name': name, 'version': version, 'kind': kind,
        'description': '测试插件',
        ('provides' if kind == 'tool' else 'databases'): mapping
        or {'demo_tool': 'bin/demo_tool.bat'},
    }
    if manifest_extra:
        body.update(manifest_extra)
    manifest = json.dumps(body, ensure_ascii=False)
    with zipfile.ZipFile(path, 'w') as z:
        if layout == 'subdir':
            z.writestr(f'{name}/plugin.json', manifest)
            z.writestr(f'{name}/bin/demo_tool.bat', '@echo demo\r\n')
            z.writestr(f'{name}/README.md', 'demo')
        else:
            z.writestr('plugin.json', manifest)
            z.writestr('bin/demo_tool.bat', '@echo demo\r\n')
        if extra_member is not None:
            z.writestr(extra_member, 'evil')
    return path


@pytest.fixture()
def plugin_root(tmp_path):
    root = tmp_path / 'plugins'
    root.mkdir()
    return str(root)


# ------------------------------------------------------------------
# 安装
# ------------------------------------------------------------------
def test_install_subdir_layout(plugin_root, tmp_path):
    zp = _make_plugin_zip(str(tmp_path / 'p.zip'))
    mf = pm.install_zip(zp, root=plugin_root)
    assert mf['name'] == 'demo-tool'
    pdir = os.path.join(plugin_root, 'demo-tool')
    assert os.path.isfile(os.path.join(pdir, 'bin', 'demo_tool.bat'))
    tools, dbs = pm.resolve_plugins(root=plugin_root)
    assert tools.get('demo_tool', '').endswith('demo_tool.bat')
    assert dbs == {}


def test_install_flat_layout(plugin_root, tmp_path):
    zp = _make_plugin_zip(str(tmp_path / 'p.zip'), layout='flat')
    mf = pm.install_zip(zp, root=plugin_root)
    assert mf['name'] == 'demo-tool'
    assert pm.resolve_plugins(root=plugin_root)[0]['demo_tool']


def test_install_rejects_missing_provides_file(plugin_root, tmp_path):
    zp = str(tmp_path / 'p.zip')
    with zipfile.ZipFile(zp, 'w') as z:
        z.writestr('demo-tool/plugin.json', json.dumps({
            'name': 'demo-tool', 'version': '1.0', 'kind': 'tool',
            'provides': {'demo_tool': 'bin/不存在的.exe'}}))
    with pytest.raises(pm.PluginError, match='不存在'):
        pm.install_zip(zp, root=plugin_root)
    # 坏包不得污染插件目录（无残留半成品目录）
    assert os.listdir(plugin_root) == []


def test_install_rejects_bad_name(plugin_root, tmp_path):
    zp = str(tmp_path / 'p.zip')
    _make_plugin_zip(zp, name='../escape')
    # 名字本身非法（..）;zip 成员也会先被穿越检查拦下——两个口子拦哪个都算对
    with pytest.raises(pm.PluginError, match='name|穿越'):
        pm.install_zip(zp, root=plugin_root)
    assert os.listdir(plugin_root) == []


def test_install_no_manifest(plugin_root, tmp_path):
    zp = str(tmp_path / 'p.zip')
    with zipfile.ZipFile(zp, 'w') as z:
        z.writestr('some/file.txt', 'no manifest here')
    with pytest.raises(pm.PluginError, match='plugin.json'):
        pm.install_zip(zp, root=plugin_root)
    assert os.listdir(plugin_root) == []


def test_install_corrupt_zip(plugin_root, tmp_path):
    zp = str(tmp_path / 'p.zip')
    with open(zp, 'wb') as f:
        f.write(b'not a zip at all')
    with pytest.raises(pm.PluginError, match='zip'):
        pm.install_zip(zp, root=plugin_root)


# ------------------------------------------------------------------
# zip 路径穿越防护
# ------------------------------------------------------------------
@pytest.mark.parametrize('member', [
    '../evil.exe',                       # 相对穿越
    'demo-tool/../../evil.exe',          # 中段穿越
    '/abs/evil.exe',                     # 绝对路径
    'C:/evil.exe',                       # 盘符
    '..\\evil.exe',                      # 反斜杠穿越
])
def test_zip_traversal_blocked(plugin_root, tmp_path, member):
    zp = _make_plugin_zip(str(tmp_path / 'p.zip'), extra_member=member)
    with pytest.raises(pm.PluginError, match='穿越|绝对路径|盘符|NUL'):
        pm.install_zip(zp, root=plugin_root)
    # 恶意成员一个都不能落盘
    for cur, _dirs, files in os.walk(plugin_root):
        assert 'evil.exe' not in files, f'穿越文件落盘: {cur}'


# ------------------------------------------------------------------
# 升级 / 卸载
# ------------------------------------------------------------------
def test_upgrade_overwrites(plugin_root, tmp_path):
    zp1 = _make_plugin_zip(str(tmp_path / 'v1.zip'), version='1.0')
    zp2 = _make_plugin_zip(str(tmp_path / 'v2.zip'), version='2.0')
    pm.install_zip(zp1, root=plugin_root)
    mf = pm.install_zip(zp2, root=plugin_root)
    assert mf['version'] == '2.0'
    assert pm.installed_plugins(plugin_root)[0]['version'] == '2.0'
    # 无升级残骸（.old / .tmp 目录清干净）
    leftovers = [e for e in os.listdir(plugin_root) if e.startswith('.')]
    assert leftovers == []


def test_uninstall(plugin_root, tmp_path):
    zp = _make_plugin_zip(str(tmp_path / 'p.zip'))
    pm.install_zip(zp, root=plugin_root)
    pm.uninstall('demo-tool', root=plugin_root)
    assert not os.path.isdir(os.path.join(plugin_root, 'demo-tool'))
    with pytest.raises(pm.PluginError, match='未安装'):
        pm.uninstall('demo-tool', root=plugin_root)
    with pytest.raises(pm.PluginError, match='非法'):
        pm.uninstall('../x', root=plugin_root)


# ------------------------------------------------------------------
# 数据库插件
# ------------------------------------------------------------------
def test_database_plugin(plugin_root, tmp_path):
    zp = str(tmp_path / 'db.zip')
    with zipfile.ZipFile(zp, 'w') as z:
        z.writestr('demo-db/plugin.json', json.dumps({
            'name': 'demo-db', 'version': '1.0', 'kind': 'database',
            'description': '测试库', 'databases': {'virus': 'db/plant'}}))
        z.writestr('demo-db/db/plant/hash_0.k2d', 'fake')
    pm.install_zip(zp, root=plugin_root)
    tools, dbs = pm.resolve_plugins(root=plugin_root)
    assert tools == {}
    assert dbs['virus'].endswith(os.path.join('db', 'plant'))


# ------------------------------------------------------------------
# Cfg 填缺合并（跳过 Config.__init__，只挂被测属性）
# ------------------------------------------------------------------
def _stub_cfg(tools, plugin_tools, plugin_databases=None, databases=None):
    from Virus_Platform_Core.config import Config
    cfg = Config.__new__(Config)
    cfg.tools = tools
    cfg.plugin_tools = plugin_tools
    cfg.plugin_databases = plugin_databases or {}
    cfg.databases = databases or {}
    return cfg


def test_cfg_tool_fills_gap(tmp_path):
    exe = tmp_path / 'demo_tool.bat'
    exe.write_text('@echo off\r\n', encoding='utf-8')
    cfg = _stub_cfg({'demo_tool': None}, {'demo_tool': str(exe)})
    assert cfg.tool('demo_tool') == str(exe)          # 探测缺 → 插件补
    assert cfg.tool_status()['demo_tool'] == str(exe)


def test_cfg_tool_platform_priority(tmp_path):
    plat = tmp_path / 'platform_ver.bat'
    plug = tmp_path / 'plugin_ver.bat'
    for p in (plat, plug):
        p.write_text('@echo off\r\n', encoding='utf-8')
    cfg = _stub_cfg({'demo_tool': str(plat)}, {'demo_tool': str(plug)})
    assert cfg.tool('demo_tool') == str(plat)         # platform.json 优先
    assert cfg.tool_status()['demo_tool'] == str(plat)


def test_cfg_tool_missing_everywhere():
    cfg = _stub_cfg({'nope': None}, {})
    with pytest.raises(FileNotFoundError, match='插件'):
        cfg.tool('nope')


def test_cfg_refresh_plugins_databases_gap_fill(tmp_path, monkeypatch):
    root = tmp_path / 'plugins'
    root.mkdir()
    zp = str(tmp_path / 'db.zip')
    with zipfile.ZipFile(zp, 'w') as z:
        z.writestr('demo-db/plugin.json', json.dumps({
            'name': 'demo-db', 'version': '1.0', 'kind': 'database',
            'databases': {'virus': 'db/plant'}}))
        z.writestr('demo-db/db/plant/hash_0.k2d', 'fake')
    pm.install_zip(zp, root=str(root))
    monkeypatch.setattr(pm, 'PLUGIN_DIR', str(root))
    cfg = _stub_cfg({}, {}, databases={'host': ''})   # host 已有（空串≠目录）、virus 缺
    cfg.refresh_plugins()
    assert 'virus' in cfg.plugin_databases
    assert cfg.databases['virus'].endswith(os.path.join('db', 'plant'))


# ------------------------------------------------------------------
# Web API（multipart 上传 → 列表 → 卸载 全链路，Flask test client）
# ------------------------------------------------------------------
@pytest.fixture()
def client(monkeypatch, tmp_path):
    monkeypatch.setattr(pm, 'PLUGIN_DIR', str(tmp_path / 'plugins'))
    import app as app_module
    app_module.app.config['TESTING'] = True
    with app_module.app.test_client() as c:
        yield c


def test_api_install_list_uninstall_roundtrip(client, tmp_path):
    zp = _make_plugin_zip(str(tmp_path / 'p.zip'), name='api-demo')
    with open(zp, 'rb') as f:
        r = client.post('/api/plugin/install',
                        data={'file': (f, 'p.zip')},
                        content_type='multipart/form-data')
    assert r.status_code == 200, r.get_data(as_text=True)
    assert r.get_json()['installed'] == 'api-demo'

    r = client.get('/api/plugin/list')
    items = r.get_json()['items']
    assert [it['name'] for it in items] == ['api-demo']
    assert items[0]['ok'] is True

    r = client.post('/api/plugin/uninstall', json={'name': 'api-demo'})
    assert r.get_json()['ok'] is True
    r = client.get('/api/plugin/list')
    assert r.get_json()['items'] == []


def test_api_install_rejects_bad_zip(client, tmp_path):
    zp = str(tmp_path / 'bad.zip')
    with open(zp, 'wb') as f:
        f.write(b'junk')
    with open(zp, 'rb') as f:
        r = client.post('/api/plugin/install',
                        data={'file': (f, 'bad.zip')},
                        content_type='multipart/form-data')
    assert r.status_code == 400
    assert r.get_json()['ok'] is False


# ------------------------------------------------------------------
# 插件源索引（fetch_index / available_from_index / sha256 校验）
# ------------------------------------------------------------------
def _write_index(tmp_path, plugins):
    idx = tmp_path / 'index.json'
    idx.write_text(json.dumps({'schema': 1, 'plugins': plugins}),
                   encoding='utf-8')
    return idx.as_uri()


def test_fetch_index_local_uri(tmp_path):
    uri = _write_index(tmp_path, [
        {'name': 'a', 'url': 'https://x/a.zip', 'version': '1'},
        {'bogus': True},                                # 脏条目被滤掉
        {'name': 'b', 'url': 'https://x/b.zip'}])
    plugins, err = pm.fetch_index(uri, force=True)
    assert err is None
    assert [p['name'] for p in plugins] == ['a', 'b']


def test_fetch_index_bad_content(tmp_path):
    uri = _write_index(tmp_path, [])
    (tmp_path / 'index.json').write_text('not json', encoding='utf-8')
    plugins, err = pm.fetch_index(uri, force=True)
    assert plugins is None and err


def test_available_from_index_detection(tmp_path, monkeypatch):
    # 本机已有 trimal（tools 表在位）、缺 fasttree；fasttree 已装同名插件目录
    root = tmp_path / 'plugins'
    (root / 'fasttree').mkdir(parents=True)
    monkeypatch.setattr(pm, 'PLUGIN_DIR', str(root))
    uri = _write_index(tmp_path, [
        {'name': 'trimal', 'url': 'https://x/trimal.zip',
         'provides': {'trimal': 'trimal.exe'}, 'size': 1},
        {'name': 'fasttree', 'url': 'https://x/ft.zip',
         'provides': {'fasttree': 'FastTree.exe'}, 'size': 2},
        {'name': 'gblocks', 'url': 'https://x/gb.zip',
         'provides': {'gblocks': 'Gblocks.exe'}, 'size': 3},
    ])
    cfg = _stub_cfg({'trimal': r'C:\x\trimal.exe'}, {})
    items, err = pm.available_from_index(cfg, index_url=uri, force=True)
    assert err is None
    by = {i['name']: i for i in items}
    assert by['trimal']['needed'] is False              # 本机已有 trimal
    assert by['fasttree']['installed'] is True          # 同名插件目录在
    assert by['gblocks']['needed'] is True              # gblocks 键本机缺
    assert by['gblocks']['missing_keys'] == ['gblocks']
    # 排序：needed 的排前（fasttree/gblocks 的键在本 stub 里都缺），trimal 殿后
    assert {items[0]['name'], items[1]['name']} == {'fasttree', 'gblocks'}
    assert items[-1]['name'] == 'trimal'


def test_install_zip_sha256_check(plugin_root, tmp_path):
    zp = _make_plugin_zip(str(tmp_path / 'p.zip'))
    good = pm._sha256_file(zp)
    with pytest.raises(pm.PluginError, match='sha256'):
        pm.install_zip(zp, root=plugin_root, expected_sha256='0' * 64)
    assert os.listdir(plugin_root) == []                # 校验失败不落盘
    mf = pm.install_zip(zp, root=plugin_root, expected_sha256=good.upper())
    assert mf['name'] == 'demo-tool'


# ------------------------------------------------------------------
# 多地址回退（国内镜像优先 → GitHub 兜底）
# ------------------------------------------------------------------
def test_install_urls_falls_back(plugin_root, tmp_path, monkeypatch):
    """第一个地址失败（连接错误）→ 自动切第二个并成功。"""
    zp = _make_plugin_zip(str(tmp_path / 'p.zip'))
    good = 'file:///' + zp.replace(os.sep, '/')

    def fake_one(url, root=None, progress=None, expected_sha256=None):
        if url == 'http://127.0.0.1:1/dead.zip':
            raise OSError('connection refused')
        return pm.install_zip(zp, root=root, expected_sha256=expected_sha256)

    monkeypatch.setattr(pm, '_install_one_url', fake_one)
    mf = pm.install_urls(['http://127.0.0.1:1/dead.zip', good], root=plugin_root)
    assert mf['name'] == 'demo-tool'


def test_install_urls_all_fail(plugin_root):
    def always_fail(*a, **kw):
        raise OSError('boom')
    orig = pm._install_one_url
    pm._install_one_url = always_fail
    try:
        with pytest.raises(pm.PluginError, match='所有下载地址均失败'):
            pm.install_urls(['http://a/x.zip', 'http://b/x.zip'], root=plugin_root)
    finally:
        pm._install_one_url = orig


def test_install_urls_empty():
    with pytest.raises(pm.PluginError, match='未提供下载地址'):
        pm.install_urls([], root='x')


# ------------------------------------------------------------------
# install_to：硬编码路径工具还原到平台原路径
# ------------------------------------------------------------------
def _make_path_plugin_zip(path, name='path-tool', install_to='3rd/tools/path-tool'):
    with zipfile.ZipFile(path, 'w') as z:
        z.writestr(f'{name}/plugin.json', json.dumps({
            'name': name, 'version': '1.0', 'kind': 'tool',
            'provides': {'path-tool': 'tool.exe'}, 'install_to': install_to}))
        z.writestr(f'{name}/tool.exe', 'MZ')
        z.writestr(f'{name}/data/x.txt', 'x')
    return path


def test_install_to_restores_original_path(tmp_path, monkeypatch):
    """install_to 型：内容摊到平台原路径，不进 plugins/<name>/。"""
    fake_root = tmp_path / 'platform'
    fake_root.mkdir()
    proot = fake_root / 'plugins'
    proot.mkdir()
    monkeypatch.setattr(pm, 'PLATFORM_ROOT', str(fake_root))
    zp = _make_path_plugin_zip(str(tmp_path / 'p.zip'))
    mf = pm.install_zip(zp, root=str(proot))
    assert mf['install_to'] == '3rd/tools/path-tool'
    tgt = fake_root / '3rd' / 'tools' / 'path-tool'
    assert (tgt / 'tool.exe').is_file()
    assert (tgt / 'data' / 'x.txt').is_file()
    assert not (tgt / 'plugin.json').exists()      # 清单不摊到目标目录
    assert not (proot / 'path-tool').exists()      # 不落 plugins/
    # 登记表 + 扫描
    assert 'path-tool' in pm._read_registry(str(proot))
    items = pm.installed_plugins(str(proot))
    assert [i['name'] for i in items] == ['path-tool']
    assert items[0]['ok'] is True
    # 卸载清目录 + 清登记
    pm.uninstall('path-tool', root=str(proot))
    assert not tgt.exists()
    assert pm._read_registry(str(proot)) == {}


def test_install_to_rejects_escape(tmp_path):
    zp = str(tmp_path / 'p.zip')
    _make_path_plugin_zip(zp, install_to='../../etc')
    with pytest.raises(pm.PluginError, match='install_to'):
        pm.install_zip(zp, root=str(tmp_path / 'plugins'))


def test_install_to_rejects_absolute(tmp_path):
    zp = str(tmp_path / 'p.zip')
    _make_path_plugin_zip(zp, install_to='C:/Windows/System32')
    with pytest.raises(pm.PluginError, match='install_to'):
        pm.install_zip(zp, root=str(tmp_path / 'plugins'))
