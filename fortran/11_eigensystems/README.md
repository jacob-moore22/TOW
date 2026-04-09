# Chapter 11: Eigensystems

Eigenvalue and eigenvector computation for symmetric and general matrices.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `jacobi` | 11.1 | Eigenvalues and eigenvectors of a real symmetric matrix |
| `eigsrt` | 11.1 | Sort eigenvectors into order by eigenvalue |
| `tred2`  | 11.2 | Householder reduction of a real symmetric matrix to tridiagonal form |
| `tqli`   | 11.3 | Eigenvalues and eigenvectors of a symmetric tridiagonal matrix |
| `balanc` | 11.5 | Balance a nonsymmetric matrix to improve eigenvalue accuracy |
| `elmhes` | 11.5 | Reduce a general matrix to upper Hessenberg form |
| `hqr`    | 11.6 | Eigenvalues of an upper Hessenberg matrix by the QR method |

## Notes

- **Symmetric matrices**: Use `jacobi` for small matrices (simple, robust) or `tred2` + `tqli` for larger ones (Householder reduction then QL iteration).
- **General (nonsymmetric) matrices**: Use `balanc` + `elmhes` + `hqr` as a pipeline.
- `eigsrt` reorders the eigenvectors/eigenvalues from `jacobi` into ascending or descending order.
