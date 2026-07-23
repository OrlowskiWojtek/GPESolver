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
            break;
        case Type::REAL_TIME:
            return "RT";
            break;
        case Type::FULL:
            return "FS";
            break;
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

class InitializationOption : public AbstractParam {
public:
    enum class Type {
        COS,              //!< Initialize with cosine function
        GAUSS,            //!< Initialize with single Gaussian
        MULTIPLE_GAUSS,   //!< Initialize with multiple Gaussians
        SETUP_GAUSS,      //!< Initialize with Gaussians spacially setupped
        FROM_BINARY_FILE, //!< Initialize from binary file
        FROM_TEXT_FILE,   //!< Initialize from text file
    };

    std::string to_string() override {
        switch (type) {
        case Type::COS:
            return "COS";
        case Type::GAUSS:
            return "GAUSS";
        case Type::MULTIPLE_GAUSS:
            return "MULTIPLE_GAUSS";
        case Type::SETUP_GAUSS:
            return "SETUP_GAUSS";
        case Type::FROM_BINARY_FILE:
            return "BINARY_FILE";
        case Type::FROM_TEXT_FILE:
            return "TEXT_FILE";
        }
        return "";
    }

    void from_string(std::string str) override {
        if (str == "COS")
            type = Type::COS;
        else if (str == "GAUSS")
            type = Type::GAUSS;
        else if (str == "MULTIPLE_GAUSS")
            type = Type::MULTIPLE_GAUSS;
        else if (str == "SETUP_GAUSS")
            type = Type::SETUP_GAUSS;
        else if (str == "BINARY_FILE")
            type = Type::FROM_BINARY_FILE;
        else if (str == "TEXT_FILE")
            type = Type::FROM_TEXT_FILE;
        else
            throw std::runtime_error("Unknown initialization strategy: " + str);
    }

    Type type = Type::MULTIPLE_GAUSS;
};

#endif
