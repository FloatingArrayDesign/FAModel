# Anchors Library

This subpackage of FAModel contains modules for anchor capacity 
calculations. There are two levels of fidelity in these models:

- Level 1 basic models are soil-dependent capacity curves for a 
  range of anchor types based on performing curve fits to 
  published information in anchor manuals and standards. 
- Level 2 intermediate models are quantitative calculations for
  suction caissons and plate anchors that account for quantitative 
  soil properties as well as their variation with embedment depth.

This plot gives an example of the capacity curves that can be 
produced by the intermediate model (holding capacity for a suction
embedded plate anchor) as a function of surface shear strength:

![Capacities](images/SEPLA_curves_small.PNG)

### Implemented level-1 model anchor and soil types

|             | Drag emb. | Suction | VLA | SEPLA |
|-------------|-----------|---------|-----|-------|
| Soft clay   | X         | X       | X   | X     |
| Medium clay | X         | X       | X   | X     |
| Hard clay   | X         |         |     |       |
| Sand        | X         |         |     |       |

### Parameters needed for level-2 anchor capacity models

|        **Anchor type** | **Suction** | **Suction** | **VLA**  | **SEPLA** |
|------------------------|-------------|-------------|----------|-----------|
|        **Soil type**   | **Clay**    | **Sand**    | **Clay** | **Clay**  |
| **Anchor parameters**  |             |             |          |           |
|        Diameter        | x           | x           |          |           |
|        Length          | x           | x           |          |           |
|        Area            |             |             | X        | X         |
|        Thickness       | ratio       | ratio       | ratio    | ratio     |
|       Embedment depth  |             |             | X        | X         |
| **Soil parameters**    |             |             |          |           |
|        gamma           | X           | X           | X        | X         |
|        Su0             | X           |             | X        | X         |
|        k               | X           |             | X        | X         |
|        alpha           | X           |             |          |           |
|        phi             |             | X           |          |           |


These models will continue to be expanded as data sources and time permit.

## Soil Classification Parameters

The soft, medium, and hard clay soil classes are distinguished by the following parameter ranges: 
| Soil Type (Clay)  | N-Value  | Effective Sat. Unit Weight, kN/m3 | Void Ratio | Natural Water Content in Sat. State, % | Undrained Shear Strength, kN/m2 |
|:-----------------:|:--------:|:---------------------------------:|:----------:|:--------------------------------------:|:-------------------------------:|
|     Very Soft     |  0 - 2   |             5.5 - 8.5             |  0.9 - 1.4 |                30 - 50                 |            0 - 12.5             |
|       Soft        |  2 - 4   |             5.5 - 8.5             |  0.9 - 1.4 |                30 - 50                 |            12.5 - 25            |
|       Medium      |  4 - 8   |             5.5 - 8.5             |  0.9 - 1.4 |                30 - 50                 |             25 - 50             |
|       Stiff       |  8 - 15  |             8.5 - 12              |    ~ 0.6   |                20 - 30                 |            50 - 100             |
|     Very Stiff    | 15 - 30  |             8.5 - 12              |    ~ 0.6   |                20 - 30                 |            100 - 200            |
|        Hard       |   < 30   |             8.5 - 12              |    ~ 0.6   |                20 - 30                 |              > 200              |


Sand can also be classified ranging from soft to hard, however only a single sand class is supported at this time. In the future, sand classes will follow the parameter ranges in the following table:

| Soil Type (sand) |  N-Value | Eff. Sat. Unit Weight, kN/m3 | Void Ratio | Natural Water Content in Sat. State, % | Eff. friction Angle | Relative density (%) |
|:----------------:|:--------:|:----------------------------:|:----------:|:--------------------------------------:|:-------------------:|:--------------------:|
|   Very   Loose   |    > 4   |            7 - 9.5           |    ~ 0.8   |                 25 - 30                |        < 30         |         < 15         |
|       Loose      |  4 - 10  |            7 - 9.5           |    ~ 0.8   |                 25 - 30                |       30 - 35       |        15 - 35       |
|     Compact      | 10 - 30  |          9.5 - 11.5          |   ~ 0.45   |                 12 - 16                |       35 - 40       |        35 - 65       |
|      Dense       | 30 - 50  |          9.5 - 11.5          |   ~ 0.45   |                 12 - 16                |       40 - 45       |       65 - 85        |
|    Very Dense    |   < 50   |          9.5 - 11.5          |   ~ 0.45   |                 12 - 16                |        > 45         |         > 85         |

Conversion functions are under development to switch between categories (level 1 anchor models) 
and quantitative soil metrics (level 2 anchor models).