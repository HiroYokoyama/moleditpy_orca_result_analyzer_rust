# Python版との同期範囲

同期元: [ORCA Result Analyzer 5.0.1](https://github.com/HiroYokoyama/moleditpy_orca_result_analyzer_plugin/tree/ee2ff6ec6ee7bc7557226b0d9f52049457101f6f)

コミット: `ee2ff6ec6ee7bc7557226b0d9f52049457101f6f`。
「同一」は改行コードと `orca_result_analyzer` → `orca_result_analyzer_rust`
のパッケージ名置換を除いた全文一致を意味します。UIや解析の変更は原則として
Python版で行い、共通ファイルをコピーして維持します。

## プラグインのPythonファイル

| ファイル（orca_result_analyzer_rust/） | 区分 | 差分・扱い |
| --- | --- | --- |
| `__init__.py` | 変更あり | メタデータの対応Python範囲のみ `<3.14`。GUI・登録処理は同一。 |
| `bond_analysis.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `charge_analysis.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `dipole_analysis.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `energy_analysis.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `energy_diag.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `force_analysis.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `freq_analysis.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `gif_export.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `gui.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `loading.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `logger.py` | Rust版専用 | 既存のRust版ローカルモジュール。共通コードから参照されない場合も保持。 |
| `mo_analysis.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `mo_compare.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `mo_energy_diag.py` | Rust版専用 | 既存のRust版ローカルモジュール。共通コードから参照されない場合も保持。 |
| `mo_engine.py` | 変更あり | BasisSetEngineの格子評価をRustへ委譲。CubeWriter・CalcWorker・設定読込と互換ヘルパーは共有。 |
| `nics_bridge.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `nmr_analysis.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `nmr_custom_ref_dialog.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `nmr_export.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `nmr_merge.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `nmr_plot.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `parser.py` | 変更あり | Rust一括解析を呼ぶラッパー。共通mixinの補助APIとNBO・スピン解析を使用。 |
| `parser_electronic.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `parser_properties.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `parser_spectra.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `parser_structure.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `property_analysis.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `scf_analysis.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `settings.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `spectrum_widget.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `tddft_analysis.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `thermal_analysis.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `traj_analysis.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `utils.py` | 同一 | パッケージ名の置換を除き全文一致。 |
| `vis.py` | 同一 | パッケージ名の置換を除き全文一致。 |

## Rust用ラッパーで維持する差分

- `parser.py`: 一括解析は `orca_parser_rs`。入力の最終ジョブ切り出し、状態リセット、
  共通mixinによるNBO・スピン解析、進捗・キャンセルの互換APIを維持します。
  ネイティブ解析の実行中は中断せず、前後と互換処理の境界でキャンセルを確認します。
- `mo_engine.py`: `BasisSetEngine` の数値評価は `orca_mo_rs`。
  `CubeWriter`、`CalcWorker`、`read_generation_settings`、`UnsupportedBasisError`
  はPython版からコピーします。`shells`・`start_idx` と正規化ヘルパーも公開API互換のため
  維持し、実際の格子評価はRustで検証します。
- `__init__.py`: 機能バージョンはPython版と同じ5.0.1です。現行PyO3ツールチェーンに
  合わせ、Rust版の対応Python範囲は3.9–3.13です。

## テストとデータ

Python版の全73テストモジュール（2,459個の `test_*` 関数）をコピーしています。
全テストコードはパッケージ名の置換以外は同一です。補助Pythonファイルとサンプル出力・
cube・NOTICEを含む全101ファイルを保持し、Rust固有のテストも追加で残します。
`tests/README.md` は同期元の説明、`tests/TESTS.md` はRust版の実行方法です。

Python版のテストを弱めたり、Rust版に合わない数値検証を削除したりせず、互換APIと
ネイティブ実装を修正して通す方針です。PythonカバレッジにはRust内部は含まれません。

## ネイティブ側の変更

- `crates/orca_parser_rs/src/parser.rs`: MO係数の連結した数値、最終MO/基底ブロック、
  UHFのスピン、IR遷移双極子微分、TDDFTの単位・最終一重項、最終NMR/振動ブロック、
  科学表記、SCFスキャン番号、未完了サイクルのエネルギー、CHELPG/RESP見出しを同期。
- `crates/orca_mo_rs/src/mo_engine.rs`: 数値計算の実装を維持。共有の球面調和関数・
  正規化・直交性・節の角度・実MO cube比較テストを適用します。
- Rust crateとプラグインのバージョンは今回5.0.1に揃えています。

## 同期チェック

`SYNC_SOURCE.json` に同期元と共有ファイル一覧を記録しています。

```bash
python scripts/check_python_sync.py --source ../moleditpy_orca_result_analyzer_plugin
```

同期元のコミット、共通Pythonファイル、共有MOクラス・関数、全テスト・fixtureを確認します。
CIでも指定コミットを取得して同じチェックを実行します。`parser.py` と
`BasisSetEngine` のRust固有部分は、全テストと数値比較で互換性を確認します。

リリース用workflow、`build.py`、Rust crate設定、同期チェック自体はRust版専用です。
ライセンスの扱いは `THIRD_PARTY_NOTICES.md` を参照してください。
