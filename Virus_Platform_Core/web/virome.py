# -*- coding: utf-8 -*-
"""Open-Virome 公共病毒组模块（平台壳页 + 同源 iframe 内嵌 SPA）。

自 app.py 拆出。三条路由：
    GET /virome              平台壳页（顶部导航 + 内嵌应用）
    GET /virome/             SPA 本体（iframe 的 src）
    GET /virome/<path>       静态资源

目录优先级：源码运行时直接用 open-virome/frontend/build（重新构建即时
生效，无需拷贝）；打包部署回退到 webapp/static/virome 副本。
"""
import os

from flask import Blueprint, render_template, send_from_directory

from Virus_Platform_Core.config import PLATFORM_ROOT
from Virus_Platform_Core.web.state import _WWW

bp = Blueprint('virome', __name__)

_VIROME_DIST = os.path.join(_WWW, 'static', 'virome')
_FRONTEND_BUILD = os.path.join(PLATFORM_ROOT, '3rd', 'open-virome',
                               'frontend', 'build')
if os.path.isdir(os.path.join(_FRONTEND_BUILD, 'static')):
    _VIROME_DIST = _FRONTEND_BUILD


@bp.route('/virome')
def page_virome():
    """平台内模块页：顶部平台导航 + 内嵌 Open-Virome 应用。

    查询参数（filters/palmprintOnly）在壳页 URL 与内嵌应用间双向同步，
    因此 /virome?filters=... 可直接收藏、分享、或由程序构造。
    """
    return render_template('virome.html')


@bp.route('/virome/')
def page_virome_app():
    """SPA 本体（iframe 的 src）。"""
    return send_from_directory(_VIROME_DIST, 'index.html')


@bp.route('/virome/<path:filename>')
def virome_asset(filename):
    # send_from_directory 本身防目录穿越；显式 check_path 与全站下载路由
    # （pages/prescreen/tools_api/download）统一口径，审计只查一处规则
    from ..utils import check_path
    full = check_path(os.path.join(_VIROME_DIST, filename))
    # 保留子目录相对路径（static/js/... 等）；basename 会把子路径剥掉导致
    # 所有嵌套资源 404，SPA 只剩空白壳
    rel = os.path.relpath(full, _VIROME_DIST).replace(os.sep, '/')
    return send_from_directory(_VIROME_DIST, rel)
