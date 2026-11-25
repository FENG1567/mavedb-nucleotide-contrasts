# GitHub 上传前整理报告

整理对象：论文候选版本 `SRR-20261124-C1` 的代码与复现材料。

## 已纳入

- 当前 BRCA2 主分析、方向/并列值审计、原始分数到映射表的离线重建代码。
- 当前 MSH2 6-TG/mock 主分析、敏感性端点和输入下载/哈希校验脚本。
- 四张主图、Figure S1–S3 的生成代码及冻结结果源表。
- Collection 冻结层级重聚合、27 套来源对账的可选审计入口和所需内部自有方法文件。
- 四张主表 Markdown/CSV、作者自有汇总结果、数据字典、科学结果约定和验证记录。
- MIT、CC BY 4.0、CITATION.cff、第三方来源和逐文件 SHA-256 清单。

## 授权边界

- 自编代码：MIT。
- 作者自有汇总结果和新生成诊断：CC BY 4.0。
- 两套 BRCA2 MaveDB 主分数表按各自记录的 CC0 条款保留来源归属。
- Ensembl CDS、GEO MSH2、BARD1、OpenSplice/COMPASS、ClinVar 等不被统一改授
  MIT 或 CC BY；详见 `THIRD_PARTY.md`、`LICENSE-DATA.md` 和
  `data/external/sources.csv`。

## 有意未做的事情

- 公开仓库：https://github.com/FENG1567/mavedb-nucleotide-contrasts；固定发布版本：v1.0.0。未分配DOI或Zenodo记录。
- 正文 Word、PDF、内部审计/模拟审稿回复、WPS 日志和服务器运行日志未放入仓库。
- MSH2 原始 GEO 文件只作为本地下载缓存存在，并被 `.gitignore` 排除；使用
  `python scripts/fetch_inputs.py --dataset MSH2_GSE162130` 获取并核验。
- BARD1/VHL 原始计数的完整重跑不是 Figure S2 提交所必需。BARD1 已启动并完成，
  VHL 因服务器临时补充 ZIP 被清理而未完成；官方 URL、bundle/member hash 和可选
  脚本均已记录。不要把该可选过程写成已完成的端到端验证。
- OpenSplice master 下载得到 HTTP 403 页面，因此没有把它当作有效数据；当前剪接
  材料仅保留固定配对和作者计算的覆盖/差异诊断。

## 上传前最后一步

在准备公开仓库的干净克隆中运行：

```bash
python tests/verify_manifest.py
python tests/verify_frozen_results.py
python tests/verify_recomputed_core.py
```

随后再次核对 `THIRD_PARTY.md` 中每个上游 URL、版本、归属和许可状态。本次更新后，
SHA-256清单覆盖实际发布的文件，排除单独下载的原始数据缓存；请在干净克隆中运行清单检查。两套
BRCA2 原始 CSV 在候选包中保持原始字节不变。
