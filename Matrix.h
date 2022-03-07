#ifndef MIRNA_MATRIX_H
#define MIRNA_MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;

class Matrix {

private:

    vector<vector<float>> mat;
    size_t rows{};
    size_t cols{};
    float init_value{};

public:
    Matrix() = default;

    /**
     * constructor for scoring matrix (uses values of float)
     * @param r
     * @param c
     * @param init_value
     */
    Matrix(size_t r, size_t c, float init_value){
        rows = r;
        cols = c;
        this->init_value = init_value;
        mat.resize(rows);
        for (auto & i : mat){
            i.resize(cols, this->init_value);
        }
    }

    /**
     * @return the value from the field inside the matrix
     */
    float& operator()(const size_t & r, const size_t & c){
        return this->mat[r][c];
    }

    /**
     * fills in the given field of the matrix
     * @param value
     */
    void fill (const size_t & r, const size_t & c, float value) {
        this->mat[r][c] = value;
    }

    /**
     * @return amount of rows
     */
    size_t get_rows() const{
        return this->rows;
    }

    /**
     * @return amount of columns
     */
    size_t get_cols() const {
        return this->cols;
    }

    void increment(const size_t & r, const size_t & c) {
        ++(this->mat[r][c]);
    }

    /**
     * @return sum of the values in one specific row
     * */
    float sumRow(const size_t & r) {
        float sum = 0;
        for(size_t c = 0; c < cols; ++c) {
            sum += this->mat[r][c];
        }
        return sum;
    }

    /**
     * @return sum of the values in one specific column
     * */
    float sumCol(const size_t & c) {
        float sum = 0;
        for(size_t r = 0; r < rows; ++r) {
            sum += this->mat[r][c];
        }
        return sum;
    }

    /**
     * @return sum of all values in the matrix
     * */
    float sumTotal () {
        float sum = 0;
        for(size_t r = 0; r < rows; ++r) {
            for(size_t c = 0; c < cols; ++c) {
                sum += this->mat[r][c];
            }
        }
        return sum;
    }

    /**
     * prints the matrix
     */
    void printMatrix(){
        const char separator    = ' ';
        const int nameWidth     = 5;

        for(size_t i = 0; i < rows; ++i){
            for(size_t j = 0; j < cols; ++j){
                cout << left << setw(nameWidth) << setfill(separator) << this->mat[i][j];
            }
            std::cout  << "\n";
        }
        std::cout  << "\n";
    }

};

#endif //MIRNA_MATRIX_H
