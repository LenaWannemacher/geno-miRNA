#ifndef MIRNA_INVALIDCHARACTER_H
#define MIRNA_INVALIDCHARACTER_H

#include <exception>
#include <string>

class InvalidCharacter : public std::exception
{
public:
    InvalidCharacter(char c);
    const char* what() const noexcept override;

private:
    std::string msg_;
};

#endif //MIRNA_INVALIDCHARACTER_H