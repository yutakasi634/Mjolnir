+++
title  = "SteepestDescent"
weight = 3000
+++

# SteepestDescent

最急降下法のためのシミュレータです。

[system]({{<relref "/docs/reference/system">}})を一つ、[forcefield]({{<relref "/docs/reference/forcefields">}})を一つ要求します。

## 例

```toml
[simulator]
type           = "SteepestDescent"
boundary_type  = "Unlimited"
precision      = "double"
delta          = 1e-4
threshold      = 1e-4
step_limit     = 1_000_000
save_step      = 100
```

## 入力

以下のパラメータを取ります。

- `type`: 文字列型
  - シミュレータの種類を指定します。このシミュレータを使う場合、`"SteepestDescent"`です。
- `boundary_type`: 文字列型
  - 境界条件の種類を指定します。具体的な大きさは[System]({{<relref "/docs/reference/system">}})で指定します。
  - `"Unlimited"`: 境界条件を設定しません。シミュレーションボックスは無限大の大きさになります。
  - `"PeriodicCuboid"`: 直方体型の周期境界条件を指定します。
- `precision`: 文字列型
  - シミュレーションに用いる浮動小数点数型の種類を指定します。
  - `"float"`: 32bit浮動小数点数を使用します。
  - `"double"`: 64bit浮動小数点数を使用します。
- `delta`: 浮動小数点数型
  - 一ステップで勾配の大きさのどれだけの割合だけ粒子を動かすかを指定します。
- `threshold`: 浮動小数点数型
  - 一ステップで最も大きく動いた粒子の動きがこれを下回った場合、収束したとして終了します。
- `step_limit`: 整数型
  - ステップ数の上限を指定します。このステップ数に達すると、収束したかに関わらず終了します。
- `save_step`: 整数型
  - 何ステップおきに状態を出力するか指定します。
  - この値に関わらず、収束して終了した場合は最終状態が最終ステップとして出力されます。
