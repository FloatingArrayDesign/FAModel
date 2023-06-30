# anchors library

This subpackage of FAModel contains modules for anchor capacity 
calculations.


|                        | **Anchor type** | **DEA** | **Suction** | **Suction** | **VLA**  | **SEPLA** |
|------------------------|-----------------|---------|-------------|-------------|----------|-----------|
|                        | **Soil type**   |         | **Clay**    | **Sand**    | **Clay** | **Clay**  |
| **Anchor parameters**  |                 |         |             |             |          |           |
|                        | D               |         | x           | x           |          |           |
|                        | L               |         | x           | x           |          |           |
|                        | A               |         |             |             | X        | X         |
|                        | t               |         | ratio       | ratio       | ratio    | ratio     |
|                        | Hs              | X       |             |             | X        | X         |
|                        | Bita            |         |             |             | X        | X         |
|                        |                 |         |             |             |          |           |
|                        |                 |         |             |             |          |           |
| **Loading parameters** |                 |         |             |             |          |           |
|                        |                 |         |             |             |          |           |
|                        | A_angle         |         | x           | x           |          |           |
|                        | F_angle         |         | x           | x           |          |           |
|                        |                 |         |             |             |          |           |
|                        |                 |         |             |             |          |           |
| **Soil parameters**    |                 |         |             |             |          |           |
|                        | gamma           |         | X           | X           | X        | X         |
|                        | Su0             |         | X           |             | X        | X         |
|                        | k               |         | X           |             | X        | X         |
|                        | alpha           |         | X           |             |          |           |
|                        | phi             |         |             | X           |          |           |
|                        | J               |         |             | X           |          |           |