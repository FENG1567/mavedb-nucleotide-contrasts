# GitHub 复现材料说明

本目录对应论文《Different DNA Variants, Same Amino Acid Substitution: A MaveDB Analysis Including the Pancreatic Cancer Susceptibility Gene BRCA2》的复现材料，
公开仓库为 https://github.com/FENG1567/mavedb-nucleotide-contrasts，固定发布版本为 **v1.0.0**。自编代码采用 **MIT**；作者自有汇总结果采用
**CC BY 4.0**。第三方原始数据、参考序列与外部代码不在这两个授权的覆盖范围内，
请保留 `THIRD_PARTY.md` 和各来源记录。

仓库已将分析源码、图形脚本、冻结汇总输入、四张主表的 Markdown/CSV、七张图的
参考 PDF、输入指纹、依赖说明和复现检查分开保存。正文 Word、内部模拟审稿与回复、
服务器运行日志和历史废弃结果不需要随 GitHub 上传。

请先读 [复现指南](docs/REPRODUCING.md) 与 [实测范围](docs/VALIDATION.md)。
“从冻结结果重绘图形”“从两套 BRCA2 分数重算”“从 MSH2 原始计数重算”是不同的
复现层级；不能仅凭图像生成成功就宣称全部原始数据流程已经验证。

主分析 BRCA2 Φ 为 **87.51%**；87.59% 只属于类别对照敏感性分析。
MSH2 主端点为 **6-TG/mock**，不是旧版本的 selected/P0。完整结果身份与适用范围见
[科学结果约定](docs/SCIENTIFIC_CONTRACT.md)。

投稿对应版本链接：https://github.com/FENG1567/mavedb-nucleotide-contrasts/releases/tag/v1.0.0。
当前没有分配软件 DOI 或论文 DOI。文件完整性清单只覆盖实际发布的文件，
按来源许可单独获取的原始数据缓存不在该发布清单内。
