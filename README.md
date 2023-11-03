# GP13 project

Welcome to GP13! These are codes written by 2 Chinese guys in LPNHE, France.

## 03/11/2023 Work

```Python
# read the selected NPZ file
npz_file = np.load('result/du_1013_threshold_5_separation_100.npz', allow_pickle=True)

# use a list comprehension to filter out empty lists
window_x = np.array([a_list for a_list in npz_file['window_x'] if len(a_list) > 0])
```
