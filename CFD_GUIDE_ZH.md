# dampmeter × OpenFOAM：單點探針量出有效阻尼 τ

**一條 scalar 時序 → 有效弛豫時間 τ（γ = 1/τ）。Zero C++ compilation — 不改 solver、不重編譯。**
三個痛點：residual 跌到 1e-5 仲敢唔敢停機？渦街／顫振衰減點量化？網格係咪將渦「抹死」？ — 呢份指南逐個答。

---

## 1. 作用限制（先讀呢段，先小人後君子）

### 量到乜 ／ 量唔到乜

| ✅ 量到 | ❌ 量唔到 |
|---|---|
| **有效阻尼時間常數 τ**（平穩窗口內嘅有效線性弛豫，γ = 1/τ） | **機制** — τ 係 effective 量，波拖／摩擦／數值黏性嘅分解要 budget 方法 |
| **τ 嘅穩定性** — 四估計量收斂 ＋ shuffle-null 對照（「個常數真係喺度」，唔係「個 fit 好靚」） | **剛度** — 勢阱曲率 ω₀² 係另一個物理量 |
| 振盪信號嘅**包絡衰減率**（Hilbert 口徑） | **流態分類** — 層流／湍流／渦街唔係 τ 直接分得出（regime 轉換只會喺 τ 變化留低間接痕跡） |
| 殘差序列嘅**時間結構**（假收斂輔助判斷） | **場結構** — 要配場算子（velocity gradient 族） |

### 證據等級

- **已驗證** ✅：大氣極渦（2,867 日，τ≈12.3 日三法收斂）、ECG、軸承振動
- **CFD：前瞻、未實測** ⚠️ — 呢個 branch 就係為此而開：*run it, break it, tell me where*

### 取樣門檻（2,800 bootstrap 實證）

| 樣本數 | 命中率 | 判 |
|---|---|---|
| 20 | 13%（τ 系統性低估 ~2×）| ❌ 硬地板 — 工具直接拒絕出數，寧可報錯唔亂猜 |
| ~180 | 57% | ⚠️ 粗層 |
| **~610** | **81.5%** | ✅ 可靠層 |
| ~1,220 | 97% | ✅ 穩定 |

注意係 610 個**採樣點**唔係 610 個時間步 — 由你嘅 `writeInterval` 決定。

### 口徑紀律（兩條鐵規）

1. **單位**：`measure()` 回傳**時間單位**（samples × dt）；`ringdown_*`／`shuffle_null_ringdown` 回傳 **samples** — 比較前必須 ×dt 對齊。
2. **兩個 ring-down 口徑唔混**：`ringdown_conditional`（衝擊條件衰減 R(k)，驗證嗰個）／`ringdown_envelope`（Hilbert 解析包絡，振盪信號用）。`|x|` 唔係包絡 — 會將 τ 量化成半個振盪周期。

### 期望管理

一次性瞬態（如 cavity spin-up）唔係平穩弛豫信號 — 四估計量好可能話 `converged: False` 甚至拒絕出數。**呢個係身份檢查做緊嘢，唔係 bug。** 真量 τ 用唔穩定算例（渦街、顫振、啟動瞬態）。

---

## 2. 使用方法（保姆級：5 分鐘接入，以官方 cavity 算例示範）

揀 OpenFOAM 最經典嘅入門算例 — **頂蓋驅動方腔流（lid-driven cavity）**。每個 OpenFOAM 用戶電腦入面都有呢個現成資料夾，唔使由零寫網格，輸出全部係純文字檔。

### Step 1：複製官方現成算例

```bash
cp -r $FOAM_TUTORIALS/incompressibleFluid/cavity ./cavity_dampmeter
cd cavity_dampmeter
```

> 版本路徑唔同：舊版喺 `$FOAM_TUTORIALS/incompressible/icoFoam/cavity`。搵唔到就喺 `$FOAM_TUTORIALS` 度 `find . -name cavity -type d`。

### Step 2：`system/controlDict` 加裝「探針」（貼 10 行字）

拉到最底，貼落去（記錄方腔正中心嘅壓力 p 同速度 U）：

```cpp
functions
{
    probes
    {
        type            probes;
        libs            ("libsampling.so");

        // ★ 關鍵：每個 timeStep 寫一次，確保採樣間隔 Δt 均勻 ★
        writeControl    timeStep;
        writeInterval   1;

        probeLocations  ( (0.05 0.05 0.005) );  // 監測點（可放多個）
        fields          (p U);
    }
}
```

順便檢查 `controlDict` 上半部：**`adjustTimeStep` 係 `no`**（cavity 默認就係），保證 Δt 固定。

### Step 3：如常生成網格跑求解器（零編譯、零修改）

```bash
blockMesh
foamRun        # 舊版 OpenFOAM 輸入 icoFoam
```

跑完目錄下多咗 `postProcessing/probes/0/p` — 純乾淨嘅「時間 vs 壓力」時序。

### Step 4：Python 讀輸出、量 τ

**快路（推薦）**：

```bash
python examples/openfoam_probes_example.py probes postProcessing/probes/0/p
```

一條命令自動完成：跳過 `#` 註釋行 → 檢查 Δt 均勻（歪咗自動重採樣＋警告）→ 去趨勢 → 四估計量出 τ → shuffle-null 檢驗。想睇 U 就 `... probes postProcessing/probes/0/U`（向量自動取模長）。

**知其所以然版（5 行核心）**：

```python
import numpy as np, dampmeter as dm

d = np.loadtxt("postProcessing/probes/0/p", comments="#")  # 列: time, probe0, ...
t, x = d[:, 0], d[:, 1]                                     # 取 1 號監測點
x = x - np.polyval(np.polyfit(np.arange(len(x)), x, 1), np.arange(len(x)))  # 去趨勢（必做）
est = dm.measure(x, dt=t[1] - t[0])
print(f"tau = {est.tau:.4g}s  converged = {est.converged}")

null = dm.shuffle_null_ringdown(x, blocks=np.arange(len(x)) // (len(x)//10))
print(f"shuffle-null tau = {null * (t[1]-t[0]):.4g}s  # 必須遠細過 tau")
```

> 想 5 分鐘睇到乾淨出數先？跑 `... demo` — 合成信號真值 τ=0.15s，實測 0.166s、3/3 收斂、null 塌縮。

### 三個現成場景，對號入座

| 你平時怎麼設 `controlDict` | 讀取哪個輸出檔 | 調用哪個方法 | 幫你解決什麼工程問題 |
|---|---|---|---|
| `probes` 監測點 | `postProcessing/probes/0/p`（或 `U`） | `dm.measure(x, dt=Δt)` | 量啟動瞬態／擾動的平復時間 τ；對比不同網格／離散格式的耗散強弱 |
| `forceCoeffs` 升阻力 | `postProcessing/forceCoeffs/0/coefficient.dat` | `dm.ringdown_envelope()` | 專抓週期性震盪（卡門渦街、機翼顫振）的**振幅包絡衰減率** |
| Solver 計算日誌 | `log.pimpleFoam` 的 residual 序列 | `dm.measure()` | **抓假收斂**：殘差表面在跌，但序列還有時間結構 ＝ 流場未真正進入物理穩態 |

### 老司機三招

**① 對比網格精度（Grid Convergence Check）**
網格加密一倍（或 `upwind` → `linearUpwind`）再跑一次，兩邊 τ_eff 對比：

```bash
python examples/openfoam_probes_example.py grid coarse/p fine/p
```

τ 明顯更短嗰套 ＝ 數值黏性更重、把流場「訓死」咗。一張 τ_eff 對比表就係「網格無關性」嘅答案。

**② 卡門渦街／升力震盪** — `probes` 換 `forceCoeffs`，讀 Cl 列：

```bash
python examples/openfoam_probes_example.py envelope postProcessing/forceCoeffs/0/coefficient.dat 2
```

**③ 假收斂熔絲（最實用）** — 每 N 步切窗口計 τ，睇佢穩唔穩：

```bash
python examples/openfoam_probes_example.py windows postProcessing/probes/0/p 0 610 305
```

residual 到 1e-5 但 **τ 跨窗劇烈跳動（CV > 50%）→ 流場未入真正穩態，千祈唔好停機**。

### 防呆警告（Gotchas）

⚠️ **開了 `adjustTimeStep yes` 必看**：Δt 在變，時間軸係歪嘅，直接算 τ 全部失真。必須先重採樣成均勻 Δt（example script 自動做＋警告）。

⚠️ **`writeInterval` 別設太疏**：跌破 20 點硬地板會被直接拒測；~610 點先可靠。

⚠️ **τ_eff ＝ 物理黏性 ＋ 數值黏性**：多數人當陷阱，但佢係**武器** — 兩套網格／格式嘅 τ_eff 對比就係量化數值耗散嘅直接方法。

⚠️ **先去趨勢再餵**：mean offset／trend 會把 ACF plateau 頂住，曲線永遠到不了 1/e，τ 系統性高估。

⚠️ **單測法不出結論**：四估計量收斂 ＋ shuffle-null 塌縮，兩個都過才算量到常數。`converged: False` 係「還沒量到」，唔係「再跑多啲就好」。

---

## 3. 技術原理（想深究先睇呢段）

### 問題設定

給定一條 scalar 時序 x（平穩窗口、有效線性弛豫），量度授動嘅 e-fold 衰減時間 τ。核心口徑：**收斂係身份檢查，唔係擬合優度** — 「呢個常數真係喺度」由獨立估計量互相印證，唔係由單一 fit 嘅 R² 話事。

### 四個獨立估計量

| 方法 | 原理 | 對應代碼 |
|---|---|---|
| `tau1_acf` | 正規化自相關 e-fold crossing（1/e 線性插值） | `acf_curve()` / `efold_from_curve()` |
| `tau2_dcor` | 距離相關 e-fold（非線性、平移不變；Székely et al. 2007 雙中心化距離矩陣） | `dcor()` |
| `tau3_ringdown` | 衝擊條件衰減 **R(k) = E[x(t)·x(t+k) | pulse] / E[x(t)² | pulse]**，pulse ＝ |Δx| > 2σ | `ringdown_conditional()` |
| `tau4_recovery` | 跌後恢復時間（輔助；右刪失傾向 — 唔可單獨作證據） | `recovery_tau()` |

**收斂規則**：≥3 個主要方法落喺中位數 ±30% 內 ⇒ 個體阻尼成立。不收斂就係「未量到」— 可能窗口唔係平穩弛豫機制，或者數據太短。

### Shuffle null（唔可跳過）

喺 block 內打亂時序（破壞時間結構、保留邊際分佈）再量 ring-down τ。**打亂後 τ 必須塌縮** — 唔塌代表你讀到嘅係數學必然（趨勢／均值 offset），唔係物理。呢個係防自欺嘅殺手鐧。輔助：`phase_surrogate_taus()`（相位隨機化、保功率譜）。

### Hilbert 解析包絡（振盪口徑）

`ringdown_envelope()`：帶通（可選）→ **Hilbert 解析包絡** |hilbert(x)| → 3 點平滑 → 正規化 MAD 門檻 `median|env − median(env)| × 1.4826` → 逐事件 log-linear 衰減 fit → 取中位數。

兩個已知陷阱（內建防護，但你要知）：`|band-passed|` 當包絡 ＝ τ 被量化成半周期；`median(env) × 1.4826` ＝ 錯嘅 MAD（門檻會浮到 ~5σ，乜都過唔到）。

### min-data 三級門檻

見 §1 表（2,800 bootstrap 特性化實證）。衝擊 ring-down 需要更多：≥40 個 pulse ≈ 1,500 樣本先有 65%。

### 代碼地圖（零黑箱）

`src/dampmeter/` 約 300 行、三個檔：`_core.py`（dcor／ACF／measure 收斂規則）、`_ringdown.py`（兩個 ring-down 口徑）、`_null.py`（shuffle null／phase surrogate）。公式逐行對應；完整推導見 repo 入面 `manuscript.pdf`（4 頁雙欄，Fig 1–3 全部真數據）。計算量：單核秒級、KB 級 RAM、每測毫秒級。

---

## 4. 還有乜嘢問題（未解決 — 歡迎帶 case 過嚟打）

1. **CFD 域未實測** ⚠️ — 呢個係最大嘅 open question，亦係呢個 branch 嘅存在意義。已驗證域係大氣／ECG／軸承；CFD 接口行唔行、邊度爛，要真 case 先知。
2. **一次性瞬態點算** — cavity 呢類 spin-up 係單次弛豫，唔係平穩窗口；工具會拒絕出數（正確行為），但「瞬態信號嘅優雅處理」未有答案。
3. **τ_eff 嘅物理／數值黏性分解** — grid study 只係間接量化；直接分解要 budget 方法，未做。
4. **非均勻 Δt 重採樣嘅偏差未量化** — 線性插值係穩陣做法，但引入幾多誤差未做特性化。
5. **residuals 用法係 heuristics** — 量嘅係時間結構唔係 residual level；完美指數衰減去趨勢後係白噪音，唔可以當「已收斂」嘅充分證據。
6. **T 算子（state detachment 紅燈）1D 接口未導出** — 極限環真收斂判定係前瞻提案，唔係現成插頭。
7. **閾值唔跨 case 搬** — `pulse_threshold=2.0` 喺唔同信號特性要重校；點樣自動校準未有方案。

---

`pip install dampmeter`｜source: github.com/MMDR10/dampmeter｜manuscript: Zenodo DOI [10.5281/zenodo.24338596](https://doi.org/10.5281/zenodo.24338596)

*AI 貢獻披露：本指南由研究者 Ma Ka Po 與 AI 研究助理（Deep Research, tygtDc）協作產出。CC BY 4.0。Judge the artifact, not the author — 有錯請帶證據過嚟。*
