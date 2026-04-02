#ifndef COMPLEX_PARAMS_HPP
#define COMPLEX_PARAMS_HPP

#include <string>

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
        SPEED_TEST      //!< Run speed test
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
        case Type::SPEED_TEST:
            return "ST";
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
        if (str == "ST")
            type = CalcStrategy::Type::SPEED_TEST;
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
        if (str == "GAUSS")
            type = Type::GAUSS;
        if (str == "MULTIPLE_GAUSS")
            type = Type::MULTIPLE_GAUSS;
        if (str == "SETUP_GAUSS")
            type = Type::SETUP_GAUSS;
        if (str == "BINARY_FILE")
            type = Type::FROM_BINARY_FILE;
        if (str == "TEXT_FILE")
            type = Type::FROM_TEXT_FILE;
    }

    Type type = Type::MULTIPLE_GAUSS;
};

class PotentialType : public AbstractParam {
public:
    enum class Type {
        REGULAR, //!< Simple potential without barier
        MEXICAN, //!< Mexican hat potential
        CRADLE,  //!< Number of minimas equal to number of cradles
    };

    std::string to_string() override {
        switch (type) {
        case Type::REGULAR:
            return "REGULAR";
            break;
        case Type::MEXICAN:
            return "MEXICAN";
            break;
        case Type::CRADLE:
            return "CRADLE";
            break;
        }

        return "";
    }

    void from_string(std::string str) override {
        if (str == "REGULAR")
            type = PotentialType::Type::REGULAR;
        if (str == "MEXICAN")
            type = PotentialType::Type::MEXICAN;
        if (str == "CRADLE")
            type = PotentialType::Type::CRADLE;
    }

    Type type = PotentialType::Type::REGULAR;
};


#endif
