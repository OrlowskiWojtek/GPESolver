#ifndef OUTPUT_FORMATTER_HPP
#define OUTPUT_FORMATTER_HPP

#include <cassert>
#include <iomanip>
#include <iostream>
#include <string>

class OutputFormatter {
public:
    static void printInfo(const std::string &&message) {
        std::cout << "[INFO]    " << message << std::endl;
    }
    static void printWarning(const std::string &&message) {
        std::cout << "[WARNING] " << message << std::endl;
    }
    static void printError(const std::string &&message) {
        std::cout << "[ERROR]   " << message << std::endl;
    }

    static void printBorderLine() {
        std::string border(box_width, '=');
        std::cout << border << std::endl;
    }

    template <typename... Args>
    static void printBoxedMessage(const Args&... args) {
        std::ostringstream oss;
        (oss << ... << args);
        std::string message = oss.str();
    
        assert(message.length() < box_width - 2 &&
               "Message too long to fit in the box.");

        if(message.length() % 2 == 0){
            message += " ";
        }

        std::cout << "|" 
                  << std::setw((box_width + message.length()) / 2) << std::right << message 
                  << std::setw((box_width - message.length()) / 2) << "|" << std::endl;
    }

private:
    static const int box_width = 60;
};

#endif
