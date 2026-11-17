# 上传步骤：GitHub 仓库 + Zenodo DOI（中文操作说明）

本文件说明如何把本包上传到 GitHub 仓库 **https://github.com/LXC27/GLP2R**，并通过 GitHub–Zenodo 集成获得 DOI。命令均为 Windows PowerShell / CMD 示例。

> 重要原则：本包内**没有**数据缓存（`data_cache/`）、R 包库（`_work/`）、中间 `.rds` 文件或任何第三方原始数据。`.gitignore` 已排除；本包目标大小约 6 MB（133 个文件）。GitHub 单文件限制 100 MB。

---

## 0. 当前状态

- 作者名单已按题名页填入 `README.md` 与 `CITATION.cff`（**6 位作者，两位通讯作者 Tianyi Tao 与 Haiquan Cao**，含单位与 ORCID），无需再改占位符。
- Zenodo 归档 DOI 已回填：**`10.5281/zenodo.24321800`**（概念 DOI `10.5281/zenodo.24321799`）；`README.md` 与 `CITATION.cff` 中不再有占位符。
- 仓库地址已填入 `https://github.com/LXC27/GLP2R`。

---

## 1. 本地检查（发布前必做）

```powershell
cd "<本仓库目录>"

# 1) 文件数与总大小（预期约 133 个文件、6 MB）
Get-ChildItem -Recurse -File | Measure-Object -Property Length -Sum |
  Select-Object Count, @{n='MB';e={[math]::Round($_.Sum/1MB,2)}}

# 2) 确认没有数据缓存/中间对象/压缩包（应无输出）
Get-ChildItem -Recurse -File -Include *.rds,*.RData,*.rda,*.gz,*.zip,*.tsv,*.bam |
  Select-Object FullName

# 3) 确认没有残留占位符（应无输出）
Select-String -Path README.md,RUNNING.md,CITATION.cff -Pattern 'PLACEHOLDER'
```

---

## 2. 初始化 git 并提交（已完成，留作流程记录）

```powershell
cd "<本仓库目录>"

git init -b main
git add .
git status                       # 逐项确认：scripts/results/tables + 8 个文档文件
git commit -m "Initial release v2.0.0: confounder-adjusted GLP2R CRC reanalysis (portable R 4.5+ scripts)"
```

检查被跟踪文件清单：

```powershell
git ls-files | Measure-Object -Line
git ls-files | Select-String -Pattern '\.(rds|RData|rda|gz|zip|tsv)$'   # 应无输出
```

---

## 3. 推送到 GitHub（已完成）

仓库 `LXC27/GLP2R` 已是 Public，`main` 为单一初始提交（133 个文件），远程 `origin` 已配置，日常只需 `git push`。首次配置时的命令：

```powershell
git remote add origin https://github.com/LXC27/GLP2R.git
git remote -v
git push -u origin main
```

### 3.1 认证

GitHub 已不支持密码推送。任选其一：

- **Git Credential Manager 弹窗**（推荐）：首次 `git push` 会弹浏览器/凭证框，登录 GitHub 授权即可。
- **Personal Access Token**：Settings → Developer settings → Personal access tokens → 生成一个（scope 勾选 `repo`），推送时以 token 作为密码；或直接把 token 写进 URL 一次：
  ```powershell
  git push https://<你的用户名>:<PAT>@github.com/LXC27/GLP2R.git main
  ```
  （不要提交含 token 的 URL 到历史，用完即撤。）

### 3.2 本机网络注意（schannel TLS 失败时）

若 `git push` 报 `schannel: AcquireCredentialsHandle failed: SEC_E_NO_CREDENTIALS`，改用 OpenSSL 后端：

```powershell
git -c http.sslBackend=openssl push -u origin main
```

后续所有 git 网络命令（`ls-remote`/`push`/`pull`）都可加 `-c http.sslBackend=openssl`。

推送后检查：打开 `https://github.com/LXC27/GLP2R`，确认文件数约 133、README 正常渲染、无数据文件；`CITATION.cff` 页右侧出现 "Cite this repository"。

---

## 4. 连接 Zenodo 并发布 Release 以获得 DOI（已完成）

1. 打开 https://zenodo.org ，用 GitHub 账号登录。
2. Settings → GitHub → Connect → 在 "GitHub repositories" 列表**打开** `LXC27/GLP2R` 的开关。
3. 回 GitHub 仓库 → Releases → Draft a new release → 打 tag → Publish。（本次实际使用的 tag 为 `v1.0.0`，Zenodo 归档文件即 `GLP2R-1.0.0.zip`；下一次发布请用 `v2.0.1`，与 `CITATION.cff` 的 2.0.0 对齐。）
4. 回 Zenodo "My uploads" 找到自动生成的记录，状态变为 Published 后复制**版本 DOI**（形如 `10.5281/zenodo.12345678`）。

Release 描述可填：

```
Reproducibility package for "Reduced GLP2R expression in colorectal cancer: a confounder-adjusted multi-cohort analysis of its prognostic, immune and epigenetic correlates".
Contents: portable R analysis scripts (43), results_master.csv, data_manifest.csv, sessionInfo.txt, result CSVs, pre-submission sensitivity outputs and manuscript tables.
Code: MIT. Derived results/tables: CC BY 4.0. Data cache and R library are not included (see README).
```

---

## 5. 回填 DOI（已完成）

版本 DOI `10.5281/zenodo.24321800`（概念 DOI `10.5281/zenodo.24321799`）已写入 `README.md`（顶部徽章 + 文末引用块）与 `CITATION.cff`（`doi:` 字段），由 `_submission_prep/GitHub/_work/fill_doi.py` 完成：

```powershell
python fill_doi.py --doi 10.5281/zenodo.24321800
```

Zenodo 已归档的快照（tag `v1.0.0`）内部仍带占位符，属正常：提交本次回填后打 `v2.0.1` tag 并发布 Release，即可让归档快照自含正确 DOI（Concept DOI 不变）。

同时把同一个 GitHub URL + DOI 填入稿件正文 "Data and code availability" 段。

---

## 6. 最终验证清单

> 已核验：133 个文件 / 5.93 MB、旧投稿痕迹为零、`https://doi.org/10.5281/zenodo.24321800` 可解析、归档 zip 内容与本仓库一致（除 DOI 占位符）。


- [ ] 仓库为 Public，文件数约 133，总大小约 6 MB；无 `data_cache/`、`_work/`、`.rds`、`.gz`、`.tsv`。
- [ ] `README.md` 渲染正常，6 位作者与两位通讯作者正确，DOI 徽章指向 `https://doi.org/<DOI>`。
- [ ] `CITATION.cff` 出现 "Cite this repository" 按钮（作者、版本、DOI、license 正确）。
- [ ] `https://doi.org/<版本DOI>` 打开 Zenodo 记录，状态 Published。
- [ ] 稿件 "Data and code availability" 填入了 GitHub URL 与 DOI。
- [ ] `Select-String -Path README.md,RUNNING.md,CITATION.cff -Pattern 'PLACEHOLDER'` 无输出。

---

## 7. 内容与文件的对应关系（便于核对）

- `scripts/`（43 个 `.R`）= 可移植分析代码（路径环境变量化 + R≥4.5 用户库分支；相对原工作区新增 `S00_install_r46.R` 与 `S01_presubmission_sensitivity.R`，不含 `100b_word_finish.ps1`）。
- `results/`（47 个文件）= `results_master.csv`（360 行，主分析数字唯一来源）、`data_manifest.csv`（65 条输入记录）、`sessionInfo.txt`、40 个结果 CSV、4 个预投稿敏感性输出。
- `tables/`（35 个文件）= 34 个表格 CSV + 合并 Excel 工作簿（预投稿版，8 基因签名）。
- 不含数据缓存（≈1.1 GB）、R 库（≈1.16 GB）、中间 `.rds`、图片文件。
- 复现限制与所需外部输入见 `README.md` 的 "What can and cannot be reproduced from this archive" 与 `RUNNING.md`。
