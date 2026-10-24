# svguest（内置副本）来源与改动说明

本目录是上游 **svguest** 的单文件构建，供平台「结果中心 → 结果图整理」使用。
浏览器内完成全部编辑与导出，图形数据不会上传到任何服务器。

| 项目 | 值 |
| --- | --- |
| 上游仓库 | https://github.com/JieqiangHe/svguest |
| 在线实例 | https://www.svguest.online/ · https://jieqianghe.github.io/svguest/ |
| 取件日期 | 2026-10-23 |
| 上游提交 | `9a7df8b9d4529d7e236356498a11b120dbe5020d`（2026-10-23） |
| 许可证 | GPL-3.0（全文见同目录 `LICENSE`） |
| 文件 | `app.html`（单文件应用，约 528 KB） |

## 本平台做的唯一改动

在 `app.html` 末尾、`</body>` 之前追加一段启动脚本：

- 支持 `?open=<同源 URL>`：页面加载完成后自动抓取该 URL 并把内容载入编辑器。
- 只接受**同源**地址；跨域 URL 直接忽略，避免把本地编辑器变成任意 URL 读取器。
- 支持 SVG（打开为可编辑矢量文档）与 PNG/JPG/GIF/WebP（作为图片放到画布上）。
- 失败时用编辑器自身的 `msg()` 提示，不影响空白画布正常使用。

除此之外未改动上游任何代码、样式或资源。

## 更新方式

```bash
# 1) 取上游最新单文件
curl -L https://www.svguest.online/app.html -o webapp/static/vendor/svguest/app.html

# 2) 重新追加 ?open= 补丁（补丁片段见 git 历史 / 本文件说明）
#    注意：上游更新后需重新核对 readFile() 与 msg() 仍是全局函数。
```

## 许可注意事项

svguest 以 GPL-3.0 发布。本平台将其作为**独立静态页面**内置，不把其代码链接进
平台 Python/JS 源码。若对外分发平台安装包，请保留本目录的 `LICENSE` 与来源说明，
并遵守 GPL-3.0 对源码提供的要求；如不希望承担该义务，可在打包时排除本目录，
改为在界面中链接上游在线实例。
