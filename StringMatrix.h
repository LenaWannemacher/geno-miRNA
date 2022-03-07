#include <utility>
#ifndef MIRNA_STRINGMATRIX_H
#define MIRNA_STRINGMATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>


class StringMatrix{

private:
    std::vector<std::vector<std::string>> str_mat;
    size_t rows{};
    size_t cols{};
    size_t all_hits{};
    std::string init_value;
    std::string mirna_identifier;

public:

    StringMatrix() = default;

    StringMatrix(size_t r, size_t c, std::string init_value){
        rows = r;
        cols = c;
        this->init_value = std::move(init_value);
        str_mat.resize(rows);
        for(auto & i : str_mat){
            i.resize(cols, this->init_value);
        }
    }

    /**
     * @return  The value from the field inside the matrix
     */
    std::string& operator()(const size_t & r, const size_t & c){
        return this->str_mat[r][c];
    }

    /**
     * Fills in the given field of the matrix
     * @param value
     */
    void fill (const size_t & r, const size_t & c, std::string& value) {
        this->str_mat[r][c] = value;
    }

    /**
     * @return  Amount of rows
     */
    size_t get_rows() const{
        return this->rows;
    }

    /**
     * @return  Amount of columns
     */
    size_t get_cols() const {
        return this->cols;
    }

    /**
     * @return  Sum of the values in one specific row
     * */
    std::string strRow(const size_t & r) {
        std::string sum;
        for(size_t c = 0; c < cols; ++c) {
            sum.append(this->str_mat[r][c]);
            sum.append(" ");
        }
        return sum;
    }

    /**
     * @return  Sum of the values in one specific column
     * */
    std::string sumCol(const size_t & c) {
        std::string sum;
        for(size_t r = 0; r < rows; ++r) {
            sum.append(this->str_mat[r][c]);
            sum.append(" ");
        }
        return sum;
    }


    /**
     * prints the matrix
     */
    void printMatrix(){
        for(size_t i = 0; i < rows; ++i){
            for(size_t j = 0; j < cols; ++j){
                std::cout << this->str_mat[i][j] << ",";
            }
            std::cout  << "\n";
        }
        std::cout  << "\n";
    }

    void setMatrix_mirna_id(std::string& mirna_id){
        mirna_identifier = mirna_id;
    }

    void setMatrix_allhits(size_t all){
        all_hits = all;
    }

    /**
     * Prints the matrix into the given output file
     *
     * @param output_file_path
     */
    void printMatrixToFile(const char* output_file_path, bool headline){

        ofstream output;
        output.open(output_file_path, std::ios_base::app);
        size_t dash_counter = 0;

        if(headline){
            output << this->mirna_identifier << endl << endl;
        }
        output << "Total hits: " << all_hits << endl;

        for(size_t r = 0; r < rows; ++r){
            for (size_t c = 0; c < cols; ++c){

                if (r % 2 == 0 || r == 0) {
                    str_mat[r][c] = '-';
                    output << right << setw(6) << setfill('-') << str_mat[r][c];
                    if(dash_counter == cols-1){
                        output << "-----";
                    }
                }
                else if(c == 0 || c == cols-1 || c%2 == 0){
                    str_mat[r][c] = '|';
                    output << str_mat[r][c];
                }
                else{
                    output << right << setw(12) << setfill(' ') << str_mat[r][c];
                }

                dash_counter++;

            }
            dash_counter = 0;
            output << "\n";
        }
        output << "\n";

        output.close();
    }

    /**
     * Prints the matrix to a given file
     *
     * @param output_file_path
     */
    void print_statmat_toFile(const char* output_file_path){

        ofstream output;
        output.open(output_file_path, std::ios_base::app);

        size_t dash_counter = 0;

        for(size_t r = 0; r < rows; ++r){
            for (size_t c = 0; c < cols; ++c){

                if (r % 2 == 0 || r == 0) {
                    str_mat[r][c] = '-';
                    output << right << setw(5) << setfill('-') << str_mat[r][c];
                    if(dash_counter == cols-1){
                        output << "-----";
                    }
                }
                else if(c == 0 || c == cols-1 || c%2 == 0){
                    str_mat[r][c] = '|';
                    output << str_mat[r][c];
                }
                else{
                    output << right << setw(12) << setfill(' ') << str_mat[r][c];
                }

                dash_counter++;

            }
            dash_counter = 0;
            output << "\n";
        }

    }

};


#endif //MIRNA_STRINGMATRIX_H
