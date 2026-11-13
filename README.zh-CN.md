# ProteinViz 中文说明

将已有的 PDB/mmCIF 蛋白模型绘制成 cartoon＋半透明表面图，支持自定义残基区间配色和旁置标签，并保存可编辑 PyMOL 会话及 JSON 参数记录。

## Windows 快速运行

建议 Python 3.11。在项目目录中运行：

```powershell
py -3.11 -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe protein_viz.py examples/data/1ubq.cif --domains @examples/domains_1ubq.json -o outputs/1ubq.png
```

安装了 uv 时也可用 `uv venv --python 3.11 .venv` 创建环境。若 PyMOL 的 pip 包不支持当前平台，请按英文 README 使用 conda-forge 环境。

## 相比早期版本修复了什么

- 修复导致程序无法启动的字符串语法错误，以及无输入时不能列出配色的问题。
- 用标准解析器读取 PDB/mmCIF，使用第一模型及作者链名、残基编号，避免依赖 mmCIF 列顺序。
- 不再把 PNG 当作结构文件读取，也不再用猜测的旋转和缩放定位标签；标注来自实际 PyMOL 视角。
- 单链自动选择，多链必须指定 `--chain`，不再硬编码 A 链。
- AgEcR1 的 603 aa 结构域范围通过 `--preset agecr1` 显式启用，不自动用于其他蛋白。
- 支持自定义区间名称；检查重叠、端点缺失和参数范围。内部缺失编号在 JSON 中记录。
- 旋转后重新适配镜头，并预留表面边界与标签区域，降低裁切和遮挡风险。
- 保存 PNG、原始 PNG、PSE、JSON；JSON 包含输入 SHA-256、链、残基、配色、视角、尺寸和 DPI。
- 移除固定 `/tmp` 路径及依赖 `which` 的查找方式，通过 PyMOL Python API 在独立实例中运行。

## 参数示例

```powershell
.venv\Scripts\python.exe protein_viz.py --list-schemes
.venv\Scripts\python.exe protein_viz.py examples/data/1ubq.cif --inspect
.venv\Scripts\python.exe protein_viz.py model.cif --chain B --domains @my_domains.json -o outputs/model.png
```

`my_domains.json` 为名称到残基起止位置的映射，例如 `{"region1": [1, 35], "region2": [36, 76]}`。提供的 ubiquitin 示例分区仅演示软件功能，不是结构域注释结论。

默认渲染已有模型，不做结构预测、结构域推断或模型置信度验证。复杂标签布局仍需人工检查。默认 300 DPI；标注图会在原始渲染宽度之外增加标签栏。API 与参数边界详见英文 README。

## 用作作品或任务素材

建议一起提供源代码、允许分享的结构文件、参数 JSON、运行环境、结果图、PSE 和复现说明。提交平台的具体要求、公开代码是否允许、AI 辅助边界和权利约定，需要按当前项目要求核对。

本次代码与文档修改有 AI 辅助，不能将它们描述为完全独立手写。此仓库原本公开并采用 MIT 许可；公开代码不应被称为未公开的独占成果。
