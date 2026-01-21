# parallel-matrix-vector-mpi
MPI implementation and performance comparison of dense and sparse (CSR) matrix-vector multiplication
# Parallel Dense vs CSR Matrix-Vector Multiplication (MPI)

## Περιγραφή
Το πρόγραμμα υλοποιεί τον πολλαπλασιασμό πίνακα–διανύσματος χρησιμοποιώντας:
- **Κανονική (dense) αναπαράσταση πίνακα**
- **Sparse αναπαράσταση CSR (Compressed Sparse Row)**

Η υλοποίηση είναι **παράλληλη** με χρήση **MPI**, και συγκρίνει:
- τον χρόνο εκτέλεσης
- τη συμπεριφορά κλιμάκωσης
- την απόδοση dense vs CSR για διαφορετικά ποσοστά μηδενικών στοιχείων.

---

## Παράμετροι εκτέλεσης
Το πρόγραμμα δέχεται 3 ορίσματα:

```bash
mpirun -np <processes> ./main <n> <percentage_of_zeros> <mult_count>
