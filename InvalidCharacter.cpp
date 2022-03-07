#include "InvalidCharacter.h"

InvalidCharacter::InvalidCharacter(char c)
        : msg_(std::string("The character ") + c + " is not valid")
{
}

const char* InvalidCharacter::what() const noexcept { return msg_.c_str(); }