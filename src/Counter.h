//
// Created by Julian Bayardo on 5/11/15.
//

#ifndef METNUM_TP2_COUNTER_H
#define METNUM_TP2_COUNTER_H

#include <string>
#include <list>
#include <map>
#include <chrono>
#include <fstream>

class Counter;
class Timer;

class Logger {
    friend class Counter;
    friend class Timer;
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
        std::cerr << name << ": " << x << std::endl;
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
        std::cerr << "RESET " << name << ": " << x << std::endl;

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
    Counter(std::string name, long long i = 0) : name(name), i(i) { }

    std::string inline getName() const {
        return this->name;
    }

    operator long long() const {
        return this->i;
    }

    Counter &operator+=(const Counter &m) {
        this->i += m.i;
        return *this;
    }

    Counter &operator-=(const Counter &m) {
        this->i -= m.i;
        return *this;
    }

    Counter &operator++() {
        this->i++;
        return *this;
    }

    Counter &operator--() {
        this->i++;
        return *this;
    }

    void set(long long x) {
        i = x;
    }

    ~Counter() {
        Logger::getInstance().reset(this->name, i);
    }
private:
    std::string name;
    long long i;
};

class Timer {
public:
    Timer(std::string name) : name(name), start(std::chrono::steady_clock::now()), end(std::chrono::steady_clock::now()), stopped(false) { }

    std::string inline getName() const {
        return this->name;
    }

    void reset(bool write = false) {
        this->end = std::chrono::steady_clock::now();

        if (write) {
            Logger::getInstance().reset(this->name, std::chrono::duration_cast<std::chrono::microseconds>(this->end - this->start).count());
        }

        this->start = std::chrono::steady_clock::now();
        this->stopped = false;
    }

    void stop() {
        if (stopped) {
            throw new std::runtime_error("Tried to stop an already stopped timer.");
        }

        this->stopped = true;
        this->end = std::chrono::steady_clock::now();
        Logger::getInstance().reset(this->name, std::chrono::duration_cast<std::chrono::microseconds>(this->end - this->start).count());
    }

    ~Timer() {
        if (!stopped) {
            this->end = std::chrono::steady_clock::now();
            Logger::getInstance().reset(this->name, std::chrono::duration_cast<std::chrono::microseconds>(this->end - this->start).count());
        }
    }
private:
    std::string name;
    std::chrono::steady_clock::time_point start;
    std::chrono::steady_clock::time_point end;
    bool stopped;
};

#endif //METNUM_TP2_COUNTER_H
