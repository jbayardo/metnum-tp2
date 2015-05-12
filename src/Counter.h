//
// Created by Julian Bayardo on 5/11/15.
//

#ifndef METNUM_TP2_COUNTER_H
#define METNUM_TP2_COUNTER_H

#include <string>
#include <list>
#include <map>
#include <fstream>

class Counter;

class Logger {
    friend class Counter;
public:
    static Logger &getInstance() {
        static Logger instance;
        return instance;
    }

    void dump(std::string file) {
        std::fstream output(file, std::ios_base::out);

        for (auto &it : this->counters) {
            output << it.first << "\t\t\t";

            for (auto &lst : it.second) {
                output << lst << " ";
            }

            output << std::endl;
        }

        output.close();
    }
private:
    void set(std::string name, long long x) {
        try {
            std::list<long long> &temporal(this->counters.at(name));

            temporal.pop_front();
            temporal.push_front(x);
        } catch(...) {
            std::list<long long> empty;
            empty.push_back(x);
            counters.insert(std::pair<std::string, std::list<long long>>(name, empty));
        }
    }

    void reset(std::string name, long long x = 0) {
        try {
            std::list<long long> &temporal(this->counters.at(name));
            temporal.push_front(x);
        } catch(...) {
            std::list<long long> empty;
            empty.push_back(x);
            counters.insert(std::pair<std::string, std::list<long long>>(name, empty));
        }
    }

    Logger() {};
    Logger(Logger const&) = delete;
    void operator=(Logger const&)  = delete;

    std::map<std::string, std::list<long long>> counters;
};

class Counter {
public:
    Counter(std::string name, long long i = 0, unsigned int writeBack = 5)
            : name(name), i(i), writeBack(writeBack), curWrite(0) {
        Logger::getInstance().reset(this->name, this->i);
    }

    std::string inline getName() const {
        return this->name;
    }

    operator long long() const {
        return this->i;
    }

    Counter &operator+=(const Counter &m) {
        this->i += m.i;
        this->curWrite++;

        if (this->curWrite >= writeBack) {
            Logger::getInstance().set(this->name, this->i);
        }

        return *this;
    }

    Counter &operator-=(const Counter &m) {
        this->i -= m.i;
        this->curWrite++;

        if (this->curWrite >= writeBack) {
            Logger::getInstance().set(this->name, this->i);
        }

        return *this;
    }

    Counter &operator++() {
        this->i++;
        this->curWrite++;

        if (this->curWrite >= writeBack) {
            Logger::getInstance().set(this->name, this->i);
        }

        return *this;
    }

    Counter &operator--() {
        this->i++;
        this->curWrite++;

        if (this->curWrite >= writeBack) {
            Logger::getInstance().set(this->name, this->i);
        }

        return *this;
    }

    void set(long long x) {
        Logger::getInstance().set(this->name, x);
    }
private:
    std::string name;
    long long i;
    unsigned int writeBack;
    unsigned int curWrite;
};

#endif //METNUM_TP2_COUNTER_H
