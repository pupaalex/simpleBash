#include "s21_matrix_oop.h"
#define S21_EPS 1e-7;

void S21Matrix::set_rows(int rows) {
  if (rows > 0 && rows != rows_) {
    S21Matrix result(rows, cols_);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols_; j++) {
        if (i >= rows_) {
          result.matrix_[i][j] = 0;
        } else {
          result.matrix_[i][j] = matrix_[i][j];
        }
      }
    }
    *this = result;
  }
}

void S21Matrix::set_cols(int cols) {
  if (cols > 0 && cols != cols_) {
    S21Matrix result(rows_, cols);
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols; j++) {
        if (j >= cols_) {
          result.matrix_[i][j] = 0;
        } else {
          result.matrix_[i][j] = matrix_[i][j];
        }
      }
    }
    *this = result;
  }
}

int S21Matrix::get_rows() { return rows_; }
int S21Matrix::get_cols() { return cols_; }

void S21Matrix::memory_allocation() {
  matrix_ = new double *[rows_]();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

// operations

bool S21Matrix::EqMatrix(const S21Matrix &other) {
  bool flag = true;
  if (rows_ == other.rows_ && cols_ == other.cols_ && matrix_ != nullptr &&
      other.matrix_ != nullptr) {
    int i = 0, j = 0;
    while (i < rows_ && flag) {
      while (j < cols_) {
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-7) {
          flag = false;
          break;
        }
        j++;
      }
      i++;
    }
  } else {
    flag = false;
  }
  return flag;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range("Error: Matrices should be the same");
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        this->matrix_[i][j] += other.matrix_[i][j];
      }
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range("Error: Matrices should be the same");
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
      }
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] * num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range("Error: Matrices should be the same");
  } else {
    S21Matrix result(cols_, other.rows_);
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < other.cols_; j++) {
        for (int k = 0; k < cols_; ++k) {
          result.matrix_[i][j] = matrix_[i][j] * other.matrix_[i][j];
        }
      }
    }
  }
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      result.matrix_[i][j] = matrix_[j][i];
    }
  }
  return result;
}

S21Matrix S21Matrix::minor_matrix(int rows, int cols) {
  S21Matrix minor(rows_ - 1, cols_ - 1);
  int row_compensator = 0;
  // корректировка индексов строк при копировании элементов
  for (int i = 0; i < minor.rows_; i++) {
    i == rows ? row_compensator = 1 : 0;
    int col_compensator = 0;
    for (int j = 0; j < minor.cols_; j++) {
      j == cols ? col_compensator = 1 : 0;
      minor.matrix_[i][j] = matrix_[i + row_compensator][j + col_compensator];
    }
  }
  return minor;
}

double S21Matrix::Determinant() {
  if (this->rows_ != this->cols_) {
    throw std::out_of_range("Error: Matrices should be the same");
  }
  double determinant = 0;
  if (this->rows_ == 1) {
    determinant = matrix_[0][0];
  } else if (this->rows_ == 2) {
    determinant = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else {
    for (int i = 0; i < rows_; i++) {
      if (matrix_[0][i] != 0) {
        double tmp = minor_matrix(0, i).Determinant();
        determinant += pow(-1, i) * matrix_[0][i] * tmp;
      }
    }
  }
  return determinant;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw std::out_of_range("Error: Matrices should be the same");
  }
  S21Matrix result(rows_, cols_);
  if (rows_ == 1 && cols_ == 1) {
    result.matrix_[0][0] = 1;
  } else {
    for (int i = 0; i < result.rows_; i++) {
      for (int j = 0; j < result.cols_; j++) {
        S21Matrix tmp = minor_matrix(i, j);
        result.matrix_[i][j] = tmp.Determinant() * pow(-1.0, i + j + 2);
      }
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = Determinant();
  if (fabs(det) < 1e-7) {
    throw std::out_of_range("Error");
  }
  S21Matrix tmp(CalcComplements().Transpose());
  tmp.MulNumber(1. / det);
  return tmp;
}

// constructors

S21Matrix::S21Matrix() {
  rows_ = 3;
  cols_ = 3;
  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows > 0 && cols > 0) {
    rows_ = rows;
    cols_ = cols;
    memory_allocation();
  } else {
    throw std::out_of_range("Invalid size of matrix");
  }
}

S21Matrix::S21Matrix(const S21Matrix &other) {
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;
  memory_allocation();
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix &&other) {
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() {
  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
  matrix_ = nullptr;
  rows_ = 0;
  cols_ = 0;
}

// operators

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const double number) {
  S21Matrix result(*this);
  result.MulNumber(number);
  return result;
}

bool S21Matrix::operator==(const S21Matrix &other) {
  return this->EqMatrix(other);
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  this->~S21Matrix();
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;
  memory_allocation();
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
  return *this;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  SumMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  SubMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  MulMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const double number) {
  MulNumber(number);
  return *this;
}

double &S21Matrix::operator()(int row, int col) {
  if (row >= rows_ || col >= cols_ || row < 0 || col < 0) {
    throw std::out_of_range("Error: Index out of range");
  } else {
    return matrix_[row][col];
  }
}