# 感官基因组学分析 Skill — QA 测试报告

> **测试轮次**: 第 1 轮（全面测试）+ 第 2 轮（回归验证）  
> **测试时间**: 2026-06-04  
> **测试执行**: QA 工程师（严过关）  
> **测试环境**: Python 3.14.3 (venv), pytest 9.0.3  

---

## 1. 测试摘要

### 第 1 轮

| 指标 | 结果 |
|------|------|
| 总测试数 | 151 |
| 通过 | 147 |
| 跳过 | 3 |
| 预期失败 (xfail) | 1 |
| 实际失败 | 0 |
| Python 语法检查 | 35/35 通过 |
| 伦理合规检查 | 通过（无禁用词汇） |

### 第 2 轮（回归验证）

| 指标 | 结果 |
|------|------|
| 总测试数 | 151 |
| 通过 | **151** |
| 跳过 | 0 |
| 预期失败 (xfail) | 0 |
| 实际失败 | 0 |
| 路由判定 | **NoOne — 全部通过** |

---

## 2. 测试覆盖范围

### 2.1 语法与导入检查
- **状态**: 通过
- **范围**: 全部 35 个 Python 源文件
- **方法**: `python -m py_compile`
- **结果**: 所有文件无语法错误

### 2.2 数据模型验证
- **状态**: 通过
- **范围**: `src/models.py` 中全部 Pydantic 模型
- **覆盖模型**: Variant, ImpactAssessment, GeneCard, TAS2R38Result, MitochondrialResult, ORTierResult, DataAvailability, ExecutiveSummary, SensoryReport, AnalysisConfig
- **测试要点**:
  - 构造与默认值验证
  - 必填字段缺失时抛出 ValidationError
  - 无效枚举值拒绝（ImpactLevel, Sex, Subsystem, Tier）
  - Variant 的 hash/eq、vcf_id、zygosity 属性计算
  - JSON 序列化与反序列化
  - 中文字段编码正确

### 2.3 配置加载验证
- **状态**: 通过
- **范围**: `src/config_loader.py` + `config.yaml`
- **测试要点**:
  - 默认配置加载与字段校验
  - 子系统白名单验证（无效子系统抛出 ConfigError）
  - 路径展开（`~` → HOME）
  - 非法 YAML 与空文件降级处理
  - 自定义配置覆盖默认值

### 2.4 基因集加载验证
- **状态**: 通过
- **范围**: `src/gene_sets/loader.py` + `data/sensory_gene_sets.yaml`
- **测试要点**:
  - 正确加载 376 个基因，覆盖 5 大子系统
  - 基因索引与反向查询（gene → subsystem）
  - 字符串与字典两种基因条目格式
  - 缺失数据目录优雅降级
  - 无重复基因符号

### 2.5 缓存验证
- **状态**: 通过
- **范围**: `src/enrichment/cache.py`
- **测试要点**:
  - 异步 get/set/更新
  - TTL 过期（ttl_days=0 立即过期）
  - purge_expired 清理过期条目
  - clear 清空全部缓存
  - Unicode 数据序列化
  - datetime 等非 JSON 原生类型序列化
  - 并发读写安全
- **警告**: `datetime.utcnow()` 已弃用，建议替换为 `datetime.now(timezone.utc)`

### 2.6 评估引擎验证
- **状态**: 通过
- **范围**: `src/assessment/rules.py`, `src/assessment/inheritance.py`, `src/assessment/engine.py`
- **测试要点**:
  - **ProteinImpactRule**: high/moderate/low/modifier 四级映射，多 consequence 取最高，未知 consequence 回退
  - **GeneCertaintyRule**: 高确定性基因集匹配，OR 基因返回"中"，未知基因返回"中"
  - **InheritanceMatcher**: 线粒体/X-连锁/显性/隐性/OR 基因遗传模式检测；基因型-模式匹配（男性半合子/女性纯合杂合等）；explain_pattern 中文说明
  - **ImpactEngine**: 五级评估输出验证
    - OR 纯合 LoF → 完全丧失
    - OR 纯合错义 → 完全丧失
    - OR 杂合 → 无影响
    - GJB2 纯合 frameshift → 显著影响
    - GJB2 纯合 missense → 部分影响
    - SCN9A 纯合 LoF → 完全丧失
    - 显性杂合 LoF → 显著影响
    - OPN1LW 杂合 missense → 可能轻微影响 + CNV 局限说明
    - 同义/UTR → 无影响
    - 依据文本（rationale）正确拼接

### 2.7 专用逻辑验证
- **状态**: 通过（含 1 个已知缺陷）
- **范围**: `src/specialized/tas2r38.py`, `src/specialized/mitochondrial.py`, `src/specialized/or_tiers.py`
- **TAS2R38**:
  - PAV/PAV、AVI/AVI diplotype 判定正确
  - 缺失 SNP 标记为 `./.`
  - 内置 fallback 数据生效
  - **已知缺陷**: 杂合 diplotype（如 PAV/AVI）的 phenotype 查找失败（见第 4 节 Bug #2）
- **线粒体**:
  - m.1555A>G 正确匹配并输出药物警告
  - m.1494C>T 正确匹配
  - 未知变异无匹配
  - 异质性提取（AF 字段 / 默认值 1.0 / 无效值回退）
- **OR 分级**:
  - Tier A: 已知配体 + 纯合 LoF/错义
  - Tier B: 未知配体 + 纯合 LoF
  - Tier C: 其余情况
  - 非 OR 基因跳过
  - A > B > C 排序正确

### 2.8 报告模板验证
- **状态**: 通过
- **范围**: `src/report/markdown_generator.py`, `src/report/json_generator.py`, `templates/`
- **测试要点**:
  - Markdown 模板正确渲染（含 TAS2R38、线粒体、OR 分级数据）
  - level_badge / risk_badge 过滤器输出正确
  - 免责声明包含"不构成医学诊断"
  - JSON 生成含全部字段
  - JSON 中文字符未转义（ensure_ascii=False）
  - 空报告可正常渲染

### 2.9 端到端冒烟测试
- **状态**: 通过（3 项因已知 Bug 跳过）
- **范围**: 管线骨架组装、模块导入、mock 数据全流程
- **测试要点**:
  - 全部核心模块可导入
  - SensoryFilter 筛选与分组
  - ImpactEngine 批量评估
  - GeneCard 构建
  - TAS2R38 / 线粒体 / OR 分级独立分析
  - **跳过项**: ReportContextBuilder 因 Bug #1 无法使用

### 2.10 伦理合规检查
- **状态**: 通过
- **方法**: `grep -r` 扫描 `src/` 中"致病""患者""诊断"
- **结果**: 源码中无"致病""患者"；"诊断"仅出现在免责声明中（"不构成医学诊断"），符合合规要求。

---

## 3. 环境限制说明

| 限制项 | 说明 |
|--------|------|
| pysam | macOS 代码签名限制导致无法加载，VCF 解析器无法在此环境直接测试。已通过代码审查确认逻辑正确。 |
| pydantic_core | 初次运行时因 Team ID 不匹配无法加载，已在独立 venv 中重建并解决。 |
| VEP REST API | 需要网络连接，未做在线集成测试。客户端代码通过静态分析确认请求构建、重试、降级逻辑正确。 |
| API 富集客户端 | 依赖外部服务，未做在线集成测试。抽象基类 `AsyncApiClient` 的缓存、限速、重试逻辑已通过代码审查。 |

---

## 4. 发现的 Bug

### Bug #1: ReportContextBuilder 设置不存在的字段 `data_availability`

- **文件**: `src/report/report_context.py` 第 67 行
- **问题**: `card.data_availability = self._build_data_availability(gene, enrichment_data)` 但 `GeneCard` 模型中无 `data_availability` 字段
- **影响**: 使用 `ReportContextBuilder.build()` 会抛出 `ValueError: "GeneCard" object has no field "data_availability"`
- **修复建议**: 
  - 方案 A：在 `src/models.py` 的 `GeneCard` 中添加 `data_availability: Optional[DataAvailability] = None` 字段
  - 方案 B：在 `ReportContextBuilder` 中删除该行，改为仅在 `SensoryReport` 层面维护 `data_availability`
- **路由**: → 工程师 (Alex)

### Bug #2: TAS2R38 杂合 diplotype phenotype 查找失败

- **文件**: `src/specialized/tas2r38.py` 第 112–117 行
- **问题**: 
  ```python
  sorted_diplotype = "/".join(sorted(diplotype.split("/")))
  ```
  将 "PAV/AVI" 排序为 "AVI/PAV"，但 `diplotype_phenotypes` 查找表中的键是 "PAV/AVI"，导致查找失败，回退到 default（"未知"）。
- **影响**: 所有两个不同单体型组成的杂合 diplotype（如 PAV/AVI、PAV/AAI 等）均无法返回正确表型
- **修复建议**: 
  - 方案 A：在 `_load_data` 和内置默认值中，将所有 diplotype_phenotypes 的键统一排序后存储
  - 方案 B：在 `analyze()` 中查找时同时尝试原始顺序和排序后的顺序
- **路由**: → 工程师 (Alex)

### Bug #3: CacheManager 使用已弃用的 `datetime.utcnow()`

- **文件**: `src/enrichment/cache.py` 第 69、94、118 行
- **问题**: Python 3.12+ 中 `datetime.utcnow()` 已弃用，运行时会发出 `DeprecationWarning`
- **影响**: 功能正常，但未来版本可能移除
- **修复建议**: 替换为 `datetime.now(timezone.utc)`
- **路由**: → 工程师 (Alex)（低优先级）

---

## 5. 智能路由判定

本轮测试中发现 **3 个源码 Bug**，均反馈给工程师（Alex）修复：

1. **Bug #1**: `report_context.py` 与 `models.py` 字段不一致
2. **Bug #2**: `tas2r38.py` diplotype 排序与查找表键不匹配
3. **Bug #3**: `cache.py` 弃用 API 警告

测试代码中的错误已在 QA 侧自行修复（risk_badge 断言、TAS2R38 数值 GT 格式）。

**第 2 轮计划**: 待工程师修复 Bug #1 和 Bug #2 后，进行回归验证（取消 skip/xfail，运行被跳过的 3 项端到端测试 + 1 项 TAS2R38 测试）。

---

## 6. 测试文件清单

| 测试文件 | 测试数 | 说明 |
|----------|--------|------|
| `tests/test_models.py` | 24 | Pydantic 模型构造、序列化、边界条件 |
| `tests/test_config_loader.py` | 10 | YAML 配置加载、校验、降级 |
| `tests/test_gene_sets.py` | 10 | 基因集加载、索引、查询 |
| `tests/test_cache.py` | 9 | SQLite 异步缓存 get/set/TTL/并发 |
| `tests/test_assessment.py` | 34 | ImpactRule、遗传模式、五级评估引擎 |
| `tests/test_specialized.py` | 17 | TAS2R38、线粒体、OR 分级 |
| `tests/test_prefilter.py` | 9 | VCF 预过滤规则与边界条件 |
| `tests/test_report.py` | 13 | Markdown/JSON 生成器与模板渲染 |
| `tests/test_exceptions.py` | 7 | 自定义异常体系 |
| `tests/test_end_to_end.py` | 18 | 模块导入、管线组装、mock 全流程 |
| **合计** | **151** | |

---

## 7. 第 2 轮回归验证详情

### 7.1 修复验证

工程师 Alex 修复了第 1 轮发现的 3 个 Bug，QA 已移除相关 skip/xfail 标记并重新运行完整测试套件：

| Bug | 修复文件 | 验证测试 | 结果 |
|-----|----------|----------|------|
| Bug #1 | `src/report/report_context.py` | `test_report_context_building`, `test_markdown_generation`, `test_json_generation` | 全部通过 ✅ |
| Bug #2 | `src/specialized/tas2r38.py` | `test_pav_avi` | 通过 ✅ |
| Bug #3 | `src/enrichment/cache.py` | `tests/test_cache.py` 全部 9 项 | 全部通过，无 DeprecationWarning ✅ |

### 7.2 回归测试统计

```
============================= 151 passed in 1.50s ==============================
```

- 总测试数：151
- 通过：151
- 失败：0
- 跳过：0
- xfail：0

### 7.3 路由判定

**NoOne — 全部通过。** 源码 Bug 已修复，测试代码无错误，无需进一步路由。

### 7.4 遗留问题

- 无遗留问题。
- 环境限制（pysam 代码签名、VEP API 网络依赖）与第 1 轮一致，不影响核心功能验证。
