// Copyright 2022 Javier Galindo

<<<<<<< HEAD
#include "../../../modules/task_2/galindo_fox_algorithm_omp/galindo_fox_algorithm_omp.h"
#include <omp.h>

bool isEqual(double x, double y) {
=======
#include "../../../modules/task_1/galindo_fox_algorithm/galindo_fox_algorithm.h"
#include <omp.h>

bool is_equal(double x, double y) {
>>>>>>> 2ee29a423ea05ad6ca04e6dbca63203a76d981d7
    return std::fabs(x - y) < 0.001;
}

bool isEqualMatrix(Matrix A, Matrix B) {
    if ((A.size() <= 0) || (B.size() <= 0))
         throw "Size of matrix must be more than 0";
    if (A.size() != B.size())
        throw "Matrices are of different size";
    size_t size = A.size();
    for (size_t i = 0; i < size; i++)
<<<<<<< HEAD
         if (!isEqual(A[i], B[i]))
=======
         if (!is_equal(A[i], B[i]))
>>>>>>> 2ee29a423ea05ad6ca04e6dbca63203a76d981d7
              return false;
    return true;
}

bool isSizeCorrect(size_t size, size_t t_count) {
    size_t blocks_count = static_cast<size_t>(sqrt(t_count));
    if (size % blocks_count != 0) {
        return false;
    }
    if (blocks_count * blocks_count != t_count) {
        return false;
    }
    return true;
}

Matrix createRandomMatrix(size_t size) {
    if (size <= 0)
          throw "Size of matrix must be more than 0";
    Matrix result(size, 0);
    std::random_device rd;
    std::mt19937 mersenne(rd());
    std::uniform_real_distribution<> urd(-100.0, 100.0);
    for (size_t i = 0; i < size; i++) {
          result[i] = static_cast<double>(urd(mersenne));
    }
    return result;
}

<<<<<<< HEAD
Matrix sequentialMatrixMultiplication(const std::vector<double>& A,
    const std::vector<double> B, size_t Size) {
=======
Matrix sequentialMatrixMultiplication(std::vector<double> A,
    std::vector<double> B, size_t Size) {
>>>>>>> 2ee29a423ea05ad6ca04e6dbca63203a76d981d7
    if (Size <= 0) {
    throw "Block size must be more than 0";
    }
    if ((A.size() <= 0) || (B.size() <= 0)) {
          throw "Size matrix must be more than 0";
    }
    if (A.size() != B.size()) {
        throw "Matrices are of different size";
    }
    if ((A.size() != Size) || (B.size() != Size)) {
        throw "Parameters and size are different";
    }
    size_t BlockSize = static_cast<size_t>(sqrt(Size));
    Matrix result(BlockSize * BlockSize, 0);
    for (size_t i = 0; i < BlockSize; i++)
        for (size_t j = 0; j < BlockSize; j++)
            for (size_t k = 0; k < BlockSize; k++)
               result[i * BlockSize + j] += A[i * BlockSize + k]
                   * B[k * BlockSize + j];
    return result;
}

<<<<<<< HEAD
Matrix sequentialBlockMatrixMultiplication(const std::vector<double>& A,
    const std::vector<double> B, size_t Size) {
=======
Matrix sequentialBlockMatrixMultiplication(std::vector<double> A,
    std::vector<double> B, size_t Size) {
>>>>>>> 2ee29a423ea05ad6ca04e6dbca63203a76d981d7
    if (Size <= 0) {
        throw "Block size must be more than 0";
    }
    if ((A.size() <= 0) || (B.size() <= 0)) {
        throw "Size matrix must be more than 0";
    }
    if (A.size() != B.size()) {
        throw "Matrices are of different size";
    }
    if ((A.size() != Size) || (B.size() != Size)) {
        throw "Parameters and size are different";
    }
    if (static_cast<size_t>(sqrt(Size)) *
        static_cast<size_t>(sqrt(Size)) != Size) {
        throw "Size is not square";
    }
    size_t BlockSize =
        static_cast<size_t>(sqrt(Size));
    size_t BlockCount =
<<<<<<< HEAD
        static_cast<size_t>(BlockSize/static_cast<size_t>(sqrt(4))) == 0
                          ? 1 : static_cast<size_t>
        (BlockSize/static_cast<size_t>(sqrt(4)));
=======
        static_cast<size_t>(BlockSize/static_cast<size_t>(sqrt(pNum))) == 0
                          ? 1 : static_cast<size_t>
        (BlockSize/static_cast<size_t>(sqrt(pNum)));
>>>>>>> 2ee29a423ea05ad6ca04e6dbca63203a76d981d7
    Matrix result(BlockSize * BlockSize, 0);
    for (size_t i = 0; i < BlockSize;
        i += BlockCount)
        for (size_t j = 0; j < BlockSize;
            j += BlockCount)
            for (size_t k = 0; k < BlockSize;
                k += BlockCount)
                for (size_t ii = i;
                    ii < ((BlockCount + i) % BlockSize + BlockCount); ii++)
                    for (size_t jj = j;
                        jj < ((BlockCount + j) % BlockSize + BlockCount); jj++)
for (size_t kk = k;
         kk < ((BlockCount + k) % BlockSize + BlockCount); kk++)
result[ii * BlockSize + jj] += A[ii * BlockSize + kk] * B[kk * BlockSize + jj];
    return result;
}

Matrix parallelBlockMatrixMultiplication(const std::vector<double>& A,
    const std::vector<double>& B, size_t Size) {
    if (Size <= 0) {
        throw "Block size must be more than 0";
    }
    if ((A.size() <= 0) || (B.size() <= 0)) {
        throw "Size matrix must be more than 0";
    }
    if (A.size() != B.size()) {
        throw "Matrices are of different size";
    }
    if ((A.size() != Size) || (B.size() != Size)) {
        throw "Parameters and size are different";
    }
    if (static_cast<size_t>(sqrt(Size)) *
        static_cast<size_t>(sqrt(Size)) != Size) {
        throw "Size is not square";
    }
    Matrix result(Size, 0);
    if (!isSizeCorrect(static_cast<size_t>(sqrt(Size)), omp_get_num_threads()))
        return result;
    size_t stage;
    size_t cols = static_cast<size_t>(sqrt(Size));
    #pragma omp parallel private(stage) shared(A, B, result)
    {
        size_t threads_count = omp_get_num_threads();
        size_t blocks_count = static_cast<size_t>(sqrt(threads_count));
        size_t block_cols_size = cols / blocks_count;
        size_t thread_num = omp_get_thread_num();
        size_t i1 = thread_num / blocks_count, j1 = thread_num % blocks_count;
        auto A1 = A.data();
        auto B1 = B.data();
        auto C1 = result.data();
        for (stage = 0; stage < blocks_count; stage++) {
            A1 = A.data() + (i1 * cols + ((i1 + stage) % blocks_count))
                * block_cols_size;
            B1 = B.data() + (((i1 + stage) % blocks_count) * cols + j1)
                * block_cols_size;
            C1 = result.data() + (i1 * cols + j1) * block_cols_size;
            for (size_t i = 0; i < block_cols_size; i++) {
                for (size_t j = 0; j < block_cols_size; j++) {
                    double tmp = 0.0;
                    for (size_t k = 0; k < block_cols_size; k++) {
                        tmp += *(A1 + i * cols + k) * *(B1 + k * cols + j);
                    }
                    *(C1 + i * cols + j) += tmp;
                }
            }
        }
    }
    return result;
}
