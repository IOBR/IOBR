# IOBR CRAN合规与v2.2.2版本对比测试设计

**日期:** 2026-05-29
**版本:** v2.2.3
**状态:** 已批准

## 概述

本设计文档描述IOBR包CRAN合规检查、网络测试增强、以及v2.2.2版本功能一致性验证的完整方案。

## 核心变更

### 1. 并行网络连通性测试

**文件:** `R/download_data.R`
**位置:** `has_internet()` 函数 (约第74-82行)

**当前状态:** 仅测试 `https://www.google.com`

**修改方案:**
```r
has_internet <- function() {
  # 同时测试两个站点，任一成功即返回TRUE
  results <- parallel::mclapply(
    c("https://www.google.com", "https://www.baidu.com"),
    function(url_str) {
      tryCatch({
        con <- url(url_str)
        on.exit(close(con), add = TRUE)
        readLines(con, n = 1)
        TRUE
      }, error = function(e) FALSE)
    },
    mc.cores = 2
  )

  # 任一成功则返回TRUE
  any(results)
}
```

**设计决策:**
- 并行测试Google和Baidu
- 使用 `parallel::mclapply` 实现并发
- 任一站点成功即返回TRUE
- 不考虑Windows兼容性（用户明确不需要）

## CRAN合规检查

### 检查要求

必须通过的测试场景:
1. **标准检查:** `devtools::check()` → 0 errors | 0 warnings | 0 notes
2. **donttest模式:** `--run-donttest` 所有example可执行
3. **断网环境:** 网络不可用时所有函数优雅返回NULL

### 执行流程

```
devtools::document() → devtools::check() → 修复问题 → 重新检查 → 通过
```

### 关键验证点

- 网络断开时函数优雅返回NULL（不抛错）
- `tempdir()` 缓存策略不写入用户目录
- 所有example使用模拟数据（无远程依赖）

## v2.2.2版本对比测试

### 测试策略

1. 从v2.2.2版本的 `.Rd` 文件提取所有modified函数的 `@examples`
2. 生成统一测试脚本 `test_v2.2.2_comparison.R`
3. 执行对比测试，发现差异立即修复
4. 目标：v2.2.2与当前版本功能100%一致

### 测试脚本结构

```r
# test_v2.2.2_comparison.R

# 1. 加载v2.2.2版本
# 2. 执行每个example，捕获输出
# 3. 切换到当前版本
# 4. 重复执行每个example
# 5. 对比两次输出
# 6. 发现差异 → 记录 → 修复 → 重新测试
```

### 对比内容

- 输出类型一致性 (data.frame, matrix, list等)
- 数值精度对比 (tolerance 1e-10)
- 结构一致性 (列名、行数、维度)
- 关键返回值对比

### 差异处理流程

```
发现差异 → 分析原因 → 修复代码 → 重新运行对比测试 → 确认一致 → 继续
```

**注意:** 测试脚本不入 `.Rbuildignore`，差异报告不commit，目标是解决差异而非仅记录。

## Git提交策略

| 步骤 | 提交内容 | Commit Message |
|------|----------|----------------|
| 1 | Baidu并行网络测试 + 文档更新 | `feat: add parallel internet test (Google + Baidu)` |
| 2 | CRAN合规修复 (如有差异) | `fix: resolve v2.2.2 compatibility differences` |

**不提交的内容:**
- 测试脚本 `test_v2.2.2_comparison.R`
- 差异报告

## 任务清单

| # | 任务 | 输出 | Commit? |
|---|------|------|---------|
| 1 | 修改 `has_internet()` | `download_data.R` | ✅ |
| 2 | 运行 `devtools::document()` | `.Rd` 文件同步 | ✅ |
| 3 | 运行 `devtools::check()` | 0 errors/warnings/notes | - |
| 4 | 运行 `--run-donttest` | 所有example通过 | - |
| 5 | 模拟断网测试 | 优雅返回NULL | - |
| 6 | 提交验证修改 | Git commit | ✅ |
| 7 | 提取v2.2.2 examples生成测试脚本 | 测试脚本 | ❌ |
| 8 | 执行对比测试 | 发现差异 | - |
| 9 | **修复差异** | 功能一致 | ✅ (如有) |
| 10 | 最终验证 | 全测试通过 | - |

## 成功标准

1. `devtools::check()` 全部通过（有网/无网/donttest三种场景）
2. v2.2.2与当前版本所有modified函数输出一致
3. 所有差异已修复并验证

## 设计决策记录

| 决策 | 选择 | 原因 |
|------|------|------|
| 网络测试策略 | 并行测试Google+Baidu | 用户偏好，最大化连通性 |
| 测试覆盖率 | 全回归测试 (100+函数) | 用户明确要求全面覆盖 |
| 测试执行方式 | 单脚本顺序执行 | 易于管理和review |
| 差异处理 | 发现即修复 | 目标是解决差异而非记录 |
| Windows兼容 | 不考虑 | 用户明确不需要 |