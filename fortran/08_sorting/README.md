# Chapter 8: Sorting

Array sorting, ranking, indexing, and equivalence class determination.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `piksrt` | 8.1 | Sort an array by straight insertion |
| `piksr2` | 8.1 | Sort two arrays by straight insertion |
| `shell`  | 8.1 | Sort an array by Shell's method |
| `sort`   | 8.2 | Sort an array by quicksort (Heapsort) |
| `sort2`  | 8.2 | Sort two arrays by quicksort |
| `indexx` | 8.4 | Construct an index for an array |
| `sort3`  | 8.4 | Sort using an index to rearrange 3 or more arrays |
| `rank`   | 8.4 | Construct a rank table for an array |
| `eclass` | 8.6 | Determine equivalence classes from a list |
| `eclazz` | 8.6 | Determine equivalence classes from a procedure |
| `qcksrt` | 8.2 | Sort an array by quicksort (alternative implementation) |
| `mdian1` | 8.5 | Find the median of an array using full sort |
| `mdian2` | 8.5 | Find the median without completely sorting |

## Notes

- For general use, `sort` (Heapsort, O(N log N) guaranteed) is recommended.
- `piksrt`/`piksr2` are O(N^2) but efficient for small N.
- `indexx` creates an index array so the original data is not moved -- useful when multiple arrays must be rearranged together (via `sort3`).
- `mdian2` is faster than `mdian1` for finding the median since it avoids a full sort.
