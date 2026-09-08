#ifndef COMPLEX_PARAMS_HPP
#define COMPLEX_PARAMS_HPP

#include <string>
#include <stdexcept>

class AbstractParam {
public:
    virtual std::string to_string()       = 0;
    virtual void from_string(std::string) = 0;
    virtual ~AbstractParam()              = default;
};

class CalcStrategy : public AbstractParam {
public:
    enum class Type {
        IMAGINARY_TIME, //!< Run imaginary time evolution
        REAL_TIME,      //!< Run real time evolution
        FULL,           //!< Run both: imaginary time and real time evolution
    };

    std::string to_string() override {
        switch (type) {
        case Type::IMAGINARY_TIME:
            return "IT";
        case Type::REAL_TIME:
            return "RT";
        case Type::FULL:
            return "FS";
        }

        return "";
    }

    void from_string(std::string str) override {
        if (str == "IT")
            type = CalcStrategy::Type::IMAGINARY_TIME;
        if (str == "RT")
            type = CalcStrategy::Type::REAL_TIME;
        if (str == "FS")
            type = CalcStrategy::Type::FULL;
    }

    Type type = CalcStrategy::Type::FULL;
};

#endif
