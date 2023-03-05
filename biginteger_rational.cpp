#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>

const long double PI = 3.14159265358979323846264338327950288419716939937510;


class BigInteger {
private:
    std::vector<short int> num;
    int cnt_digits = 4;
    int radix = 10000;
    bool IsNegative = false;

    static bool operator_greater_for_abs(const BigInteger& biginteger1, const BigInteger& biginteger2) {
        if (biginteger1.num.size() != biginteger2.num.size())
            return (biginteger1.num.size() > biginteger2.num.size());
        bool IsGreater = false;
        bool IsEqual = true;
        for (int i = biginteger1.num.size() - 1; i >= 0; --i) {
            if (biginteger1.num[i] != biginteger2.num[i]) {
                IsEqual = false;
                if (biginteger1.num[i] > biginteger2.num[i])
                    IsGreater = true;
                break;
            }
        }
        return (IsEqual || IsGreater);
    }

    static bool equivalent(const BigInteger& biginteger1, const BigInteger& biginteger2) {
        if (biginteger1.num.size() != biginteger2.num.size())
            return false;
        for (size_t i = 0; i < biginteger1.num.size(); ++i)
            if (biginteger1.num[i] != biginteger2.num[i])
                return false;
        return true;
    }


    void fft(std::vector<std::complex<long double> >& cur, int rev, int mx) const {
        std::complex<long double> temp1[mx / 2], temp2[mx / 2];
        for (int len = mx / 2; len >= 1; len /= 2) {
            for (int i = 0; i < len; ++i){
                int fl = 0;
                int pos = 0;
                for (int j = i; j < mx; j += len) {
                    if (!fl)
                        temp1[pos] = cur[j];
                    else
                        temp2[pos] = cur[j], ++pos;
                    fl = 1 - fl;
                }
                std::complex<long double> ang(std::cos(rev * 2 * PI / (mx / len)), std::sin(rev * 2 * PI / (mx / len)));
                std::complex<long double> cur_ang(1, 0);
                for (int j = i, cnt = 0; cnt < pos; j += len, ++cnt) {
                    cur[j] = temp1[cnt] + cur_ang * temp2[cnt];
                    cur[j + mx / 2] = temp1[cnt] - cur_ang * temp2[cnt];
                    cur_ang *= ang;
                }
            }
        }
    }


public:

    friend bool operator<(const BigInteger&, const BigInteger&);
    friend BigInteger operator/(const BigInteger&, const BigInteger&);
    friend BigInteger operator*(const BigInteger&, const BigInteger&);

    BigInteger() {}

    BigInteger(int number) {
        if (number < 0) {
            IsNegative = true;
            number = -number;
        }
        if (number == 0)
            num.push_back(0);
        while (number != 0) {
            num.push_back(number % 10);
            number /= 10;
        }
    }

    BigInteger(const BigInteger& biginteger) : num(biginteger.num), IsNegative(biginteger.IsNegative) {}

    BigInteger& operator+=(const BigInteger& biginteger) {
        if (((int)biginteger.IsNegative + (int)IsNegative) % 2) {
            IsNegative ^= 1;
            *this -= biginteger;
            IsNegative ^= 1;
            return *this;
        }
        int RadixOverflow = 0;
        size_t sz = std::max(num.size(), biginteger.num.size());
        for (size_t i = 0; i < sz; ++i) {
            int new_num = RadixOverflow + (i < num.size() ? num[i] : 0) + (i < biginteger.num.size() ? biginteger.num[i] : 0);
            if (i < num.size())
                num[i] = new_num % 10;
            else
                num.push_back(new_num % 10);
            RadixOverflow = new_num / 10;
        }
        if (RadixOverflow)
            num.push_back(RadixOverflow);
        update_null();
        return *this;
    }

    BigInteger& operator-=(const BigInteger& biginteger) {
        if (((int)biginteger.IsNegative + (int)IsNegative) % 2) {
            IsNegative ^= 1;
            *this += biginteger;
            IsNegative ^= 1;
            return *this;
        }
        bool IsGreaterAbs = operator_greater_for_abs(*this, biginteger);
        while(num.size() < biginteger.num.size())
            num.push_back(0);
        bool RadixOverflow = false;
        if (!IsGreaterAbs)
            IsNegative = (1 ^ biginteger.IsNegative);
        for (size_t i = 0; i < num.size(); ++i) {
            int new_num = (!IsGreaterAbs) ? biginteger.num[i] - RadixOverflow - num[i] : num[i] - RadixOverflow - (i < biginteger.num.size() ? biginteger.num[i] : 0);
            RadixOverflow = 0;
            if (new_num < 0)
                new_num += 10, RadixOverflow = 1;
            num[i] = new_num;
        }
        while (num.size() > 1 && num.back() == 0)
            num.pop_back();
        update_null();
        return *this;
    }

    BigInteger& operator%=(const BigInteger& biginteger) {
        return (*this) -= (*this) / biginteger * biginteger;
    }

    BigInteger& operator/=(const BigInteger& a) {
        BigInteger pp(*this), pp1(a);
        pp1.IsNegative = pp.IsNegative = false;
        std::reverse(pp.num.begin(), pp.num.end());
        BigInteger ans, cur;
        ans.num.push_back(0);
        bool fl = false;
        BigInteger x(pp1);
        for (size_t i = 0; i < pp.num.size(); ++i) {
            if (cur.num.size() == 1 && cur.num[0] == 0)
                cur.num.pop_back();
            cur.num.push_back(pp.num[i]);
            std::reverse(cur.num.begin(), cur.num.end());
            if (operator_greater_for_abs(cur, pp1)) {
                fl = true;
                for (int j = 2; j <= 10; ++j){
                    x += pp1;
                    if (equivalent(x, cur)) {
                        ans.num.push_back(j);
                        cur -= x;
                        for (int k = 0; k < j - 1; ++k)
                            x -= pp1;
                        break;
                    }
                    if (j == 10 || operator_greater_for_abs(x, cur)) {
                        ans.num.push_back(j - 1);
                        cur -= (x -= pp1);
                        for (int k = 0; k < j - 2; ++k)
                            x -= pp1;
                        break;
                    }
                }
            } else {
                if (fl)
                    ans.num.push_back(0);
            }
            std::reverse(cur.num.begin(), cur.num.end());
        }
        std::reverse(ans.num.begin(), ans.num.end());
        while (ans.num.size() > 1 && ans.num.back() == 0)
            ans.num.pop_back();
        num.resize(ans.num.size());
        for (size_t i = 0; i < ans.num.size(); ++i)
            num[i] = ans.num[i];
        IsNegative ^= a.IsNegative;
        update_null();
        return *this;
    }

    BigInteger& operator*=(const BigInteger& a) {
        std::vector<std::complex<long double> > pp, pp1;
        int k = 0;
        std::complex<long double> x(0, 0);
        for (size_t i = 0; i < a.num.size(); ++i) {
            ++k;
            std::complex<long double> t(a.num[i], 0);
            int tt = 1;
            for (size_t j = 0; j < (size_t)(((k % cnt_digits) - 1 + cnt_digits) % cnt_digits); ++j)
                tt *= 10;
            std::complex<long double> xx(tt, 0);
            x += xx * t;
            if (k % cnt_digits == 0 || i == a.num.size() - 1) {
                pp.push_back(x);
                std::complex<long double> xz(0, 0);
                x = xz;
            }
        }
        k = 0;
        std::complex<long double> xz(0, 0);
        x = xz;
        for (size_t i = 0; i < num.size(); ++i) {
            ++k;
            std::complex<long double> t(num[i], 0);
            int tt = 1;
            for (size_t j = 0; j < (size_t)(((k % cnt_digits) - 1 + cnt_digits) % cnt_digits); ++j)
                tt *= 10;
            std::complex<long double> xx(tt, 0);
            x += xx * t;
            if (k % cnt_digits == 0 || i == num.size() - 1) {
                pp1.push_back(x);
                std::complex<long double> xz(0, 0);
                x = xz;
            }
        }
        int mx = std::max(pp.size(), pp1.size());
        mx *= 2;
        int t = 0;
        while((1 << t) < mx)
            ++t;
        mx = (1 << t);
        while (pp.size() != (size_t)mx) {
            std::complex<long double> x(0, 0);
            pp.push_back(x);
        }
        while (pp1.size() != (size_t)mx) {
            std::complex<long double> x(0, 0);
            pp1.push_back(x);
        }
        fft(pp, 1, mx);
        fft(pp1, 1, mx);
        for (int i = 0; i < mx; ++i)
            pp[i] *= pp1[i];
        fft(pp, -1, mx);
        std::complex<long double> xx(mx, 0);
        x = xx;
        std::vector<char> ans;
        long long fl = 0;
        for (size_t i = 0; i < (size_t)mx; ++i) {
            long long t = round((pp[i] / x).real());
            fl += t;
            std::string s = "";
            int k = fl % radix;
            while (k != 0) {
                s += char('0' + k % 10);
                k /= 10;
            }
            while (s.length() != (size_t)cnt_digits)
                s += '0';
            for (size_t j = 0; j < s.length(); ++j)
                ans.push_back(s[j]);
            fl /= radix;
        }
        while (ans.size() > 1 && ans.back() == '0')
            ans.pop_back();
        num.resize(ans.size());
        for (size_t i = 0; i < ans.size(); ++i)
            num[i] = ans[i] - '0';
        IsNegative ^= a.IsNegative;
        update_null();
        return *this;
    }

    std::string division(const BigInteger& a, int k) {
        BigInteger pp(*this), pp1(a);
        std::string s = "";
        if ((*this) && (IsNegative ^ a.IsNegative))
            s += '-';
        pp1.IsNegative = pp.IsNegative = false;
        std::reverse(pp.num.begin(), pp.num.end());
        BigInteger ans, cur;
        ans.num.push_back(0);
        bool fl = false;
        int t = 0;
        BigInteger x(pp1);
        for (size_t i = 0; i < pp.num.size() || t <= k; ++i) {
            if (cur.num.size() == 1 && cur.num[0] == 0)
                cur.num.pop_back();
            if (i < pp.num.size())
                cur.num.push_back(pp.num[i]);
            else {
                ++t;
                cur.num.push_back(0);
            }
            std::reverse(cur.num.begin(), cur.num.end());
            if (operator_greater_for_abs(cur, pp1)) {
                fl = true;
                for (int j = 2; j <= 10; ++j) {
                    x += pp1;
                    if (equivalent(x, cur)) {
                        cur -= x;
                        for (int k = 0; k < j - 1; ++k)
                            x -= pp1;
                        s += char(j + '0');
                        break;
                    }
                    if (j == 10 || operator_greater_for_abs(x, cur)) {
                        s += char(j - 1 + '0');
                        cur -= (x -= pp1);
                        for (int k = 0; k < j - 2; ++k)
                                x -= pp1;
                        break;
                    }
                }
            } else {
                s += '0';
            }
            std::reverse(cur.num.begin(), cur.num.end());
        }
        fl = 0;
        if (s.back() - '0' >= 5)
            fl = 1;
        s.pop_back();
        std::string s1 = "";
        while (s.size()) {
            if ((int)s1.size() == k && k != 0)
                s1 += '.';
            if (s.back() == '-') {
                s1 += s.back();
                s.pop_back();
                continue;
            }
            int k = s.back() - '0' + fl;
            s1 += char(k % 10 + '0');
            fl = k / 10;
            s.pop_back();
        }
        if (s1.back() == '-') {
            s += '-';
            s1.pop_back();
        }
        while (s1.size() > 1 && s1.back() != '.' && s1.back() == '0')
            s1.pop_back();
        if (s1.size() > 0 && s1.back() == '.')
            s1 += '0';
        for (int i = s1.length() - 1; i >= 0; --i)
            s += s1[i];
        return s;
    }


    BigInteger& operator++() {
        return *this += 1;
    }

    BigInteger& operator--() {
        return *this -= 1;
    }

    BigInteger operator-() const {
        BigInteger a(*this);
        a.IsNegative ^= true;
        a.update_null();
        return a;
    }

    BigInteger operator++(int) {
        BigInteger x(*this);
        (*this) += 1;
        return x;
    }

    BigInteger operator--(int) {
        BigInteger x(*this);
        (*this) -= 1;
        return x;
    }

    explicit operator bool() const {
        return !(num.size() == 1 && num[0] == 0);
    }

    explicit operator int() const {
        int x = 0;
        for (int i = num.size() - 1; i >= 0; --i)
            x = x * 10 + num[i];
        if (IsNegative)
            x *= -1;
        return x;
    }

    std::string toString() const {
        std::string ret_str = "";
        if (IsNegative)
            ret_str += '-';
        for (int i = num.size() - 1; i >= 0; --i){
            ret_str += char(num[i] + '0');
        }
        return ret_str;
    }


    friend std::ostream& operator<<(std::ostream& out, const BigInteger& biginteger) {
        if (biginteger.IsNegative)
            out << '-';
        for (int i = biginteger.num.size() - 1; i >= 0; --i)
            out << biginteger.num[i];
        return out;
    }

    friend std::istream& operator>>(std::istream& in, BigInteger& biginteger) {
        biginteger.clear();
        char symbol;
        bool is_read_space_before = false;
        while(true) {
            symbol = in.get();
            if ((is_read_space_before && isspace(symbol)) || symbol == EOF)
                break;
            if (isspace(symbol))
                continue;
            is_read_space_before = true;
            if (symbol == '-')
                biginteger.IsNegative = true;
            else
                biginteger.num.push_back(symbol - '0');
        }
        std::reverse(biginteger.num.begin(), biginteger.num.end());
        return in;
    }

    void clear() {
        num.resize(0);
        IsNegative = false;
    }

    void update_null() {
        if (num.size() == 1 && num[0] == 0)
            IsNegative = false;
    }

    void change_IsNegative(bool x){
        IsNegative ^= x;
    }

    void assign_IsNegative(bool x) {
        IsNegative = x;
    }

    bool get_IsNegative() const {
        return IsNegative;
    }

    ~BigInteger() {
        num.resize(0);
    }

};




BigInteger operator-(const BigInteger& biginteger1, const BigInteger& biginteger2) {
    BigInteger x(biginteger1);
    return x -= biginteger2;
}

BigInteger operator+(const BigInteger& biginteger1, const BigInteger& biginteger2) {
    BigInteger x(biginteger1);
    return x += biginteger2;
}


BigInteger operator*(const BigInteger& biginteger1, const BigInteger& biginteger2) {
    BigInteger x(biginteger1);
    return x *= biginteger2;
}

BigInteger operator/(const BigInteger& biginteger1, const BigInteger& biginteger2) {
    BigInteger x(biginteger1);
    return x /= biginteger2;
}

BigInteger operator%(const BigInteger& biginteger1, const BigInteger& biginteger2) {
    BigInteger x(biginteger1);
    return x %= biginteger2;
}

bool operator<(const BigInteger& biginteger1, const BigInteger& biginteger2) {
    if (biginteger1.IsNegative && !biginteger2.IsNegative)
            return true;
    if (!biginteger1.IsNegative && biginteger2.IsNegative)
        return false;
    if (biginteger1.num.size() != biginteger2.num.size())
        return (biginteger1.num.size() < biginteger2.num.size()) ^ biginteger1.IsNegative;
    bool IsLess = true;
    bool IsEqual = true;
    for (int i = biginteger1.num.size() - 1; i >= 0; --i){
        if (biginteger1.num[i] != biginteger2.num[i]){
            IsEqual = false;
            if (biginteger1.num[i] > biginteger2.num[i])
                IsLess = false;
            break;
        }

    }
    if (IsEqual)
        return false;
    return IsLess ^ biginteger1.IsNegative;
}

bool operator>(const BigInteger& biginteger1, const BigInteger& biginteger2) {
    return -biginteger1 < -biginteger2;
}

bool operator<=(const BigInteger& biginteger1, const BigInteger& biginteger2) {
    return !(biginteger2 < biginteger1);
}

bool operator>=(const BigInteger& biginteger1, const BigInteger& biginteger2) {
    return !(biginteger1 < biginteger2);
}

bool operator!=(const BigInteger& biginteger1, const BigInteger& biginteger2) {
    return (biginteger1 < biginteger2 || biginteger2 < biginteger1);
}

bool operator==(const BigInteger& biginteger1, const BigInteger& biginteger2) {
    return !(biginteger1 < biginteger2 || biginteger2 < biginteger1);
}




class Rational{
private:
    BigInteger numerator;
    BigInteger denominator;
public:

    friend bool operator<(const Rational&, const Rational&);

    Rational() {}

    Rational(int x) : numerator(x), denominator(1) {}

    Rational(const BigInteger& biginteger) : numerator(biginteger), denominator(1) {}

    Rational(const Rational& rational) : numerator(rational.numerator), denominator(rational.denominator) {}

    BigInteger gcd(BigInteger biginteger1, BigInteger biginteger2) const {
        if (!biginteger2)
            return biginteger1;
        return gcd(biginteger2, biginteger1 % biginteger2);
    }

    void reduce() {
        numerator.change_IsNegative(denominator.get_IsNegative());
        denominator.change_IsNegative(denominator.get_IsNegative());
        bool IsNeg = numerator.get_IsNegative();
        numerator.assign_IsNegative(false);
        BigInteger gcd_biginteger(gcd(denominator, numerator));
        numerator.change_IsNegative(IsNeg);
        numerator /= gcd_biginteger;
        denominator /= gcd_biginteger;
    }



    Rational& operator+=(const Rational& rational) {
        numerator.change_IsNegative(denominator.get_IsNegative());
        denominator.change_IsNegative(denominator.get_IsNegative());
        numerator *= rational.denominator;
        numerator += rational.numerator * denominator;
        denominator *= rational.denominator;
        reduce();
        return *this;
    }



    Rational& operator-=(const Rational& rational) {
        numerator.change_IsNegative(denominator.get_IsNegative());
        denominator.change_IsNegative(denominator.get_IsNegative());
        numerator *= rational.denominator;
        numerator -= rational.numerator * denominator;
        denominator *= rational.denominator;
        reduce();
        return *this;
    }


    Rational& operator*=(const Rational& rational) {
        numerator.change_IsNegative(denominator.get_IsNegative());
        denominator.change_IsNegative(denominator.get_IsNegative());
        numerator *= rational.numerator;
        denominator *= rational.denominator;
        reduce();
        return *this;
    }



    Rational& operator/=(const Rational& rational) {
        numerator.change_IsNegative(denominator.get_IsNegative());
        denominator.change_IsNegative(denominator.get_IsNegative());
        numerator *= rational.denominator;
        denominator *= rational.numerator;
        reduce();
        return *this;
    }


    Rational operator-() const {
        Rational rational(*this);
        rational.numerator.change_IsNegative(true);
        return rational;
    }

    explicit operator double() {
        std::string to_str_double = asDecimal(18);
        long double ans = 0;
        long double DecRadix = 1;
        bool IsMantissa = 0;
        for (size_t i = 0; i < to_str_double.length(); ++i) {
            if (to_str_double[i] == '.') {
                IsMantissa = true;
                continue;
            }
            if (!IsMantissa)
                ans = ans * 10 + to_str_double[i] - '0';
            else {
                DecRadix *= 10;
                ans += (to_str_double[i] - '0') * (1.0 / DecRadix);
            }
        }
        if (numerator.get_IsNegative() ^ denominator.get_IsNegative())
            ans *= -1;
        return ans;
    }

    std::string toString() {
        reduce();
        std::string ret_str = "";
        numerator.assign_IsNegative(numerator.get_IsNegative() ^ denominator.get_IsNegative());
        denominator.assign_IsNegative(false);
        ret_str += numerator.toString();
        if (!numerator || !(denominator - 1))
            return ret_str;
        ret_str += '/';
        ret_str += denominator.toString();
        return ret_str;
    }

    std::string asDecimal(int k = 0) {
        return numerator.division(denominator, k);
    }



};

Rational operator-(const Rational& rational1, const Rational& rational2) {
    Rational x(rational1);
    return x -= rational2;
}

Rational operator+(const Rational& rational1, const Rational& rational2) {
    Rational x(rational1);
    return x += rational2;
}


Rational operator*(const Rational& rational1, const Rational& rational2) {
    Rational x(rational1);
    return x *= rational2;
}


Rational operator/(const Rational& rational1, const Rational& rational2) {
    Rational x(rational1);
    return x /= rational2;
}


bool operator<(const Rational& rational1, const Rational& rational2) {
    if ((rational1.numerator.get_IsNegative() ^ rational1.denominator.get_IsNegative()) && !(rational2.numerator.get_IsNegative() ^ rational2.denominator.get_IsNegative()))
            return true;
    if (!(rational1.numerator.get_IsNegative() ^ rational1.denominator.get_IsNegative()) && (rational2.numerator.get_IsNegative() ^ rational2.denominator.get_IsNegative()))
        return false;
    BigInteger bg1;
    bg1 = rational1.numerator * rational2.denominator;
    bg1.assign_IsNegative(rational1.numerator.get_IsNegative() ^ rational1.denominator.get_IsNegative());
    if (rational2.denominator.get_IsNegative())
        bg1.change_IsNegative(true);
    BigInteger bg2;
    bg2 = rational2.numerator * rational1.denominator;
    bg2.assign_IsNegative(rational2.numerator.get_IsNegative() ^ rational1.denominator.get_IsNegative());
    if (rational1.denominator.get_IsNegative())
        bg2.change_IsNegative(true);
    return bg1 < bg2;
}

bool operator>(const Rational& rational1, const Rational& rational2) {
    return -rational1 < -rational2;
}

bool operator<=(const Rational& rational1, const Rational& rational2) {
    return !(rational2 < rational1);
}

bool operator>=(const Rational& rational1, const Rational& rational2) {
    return !(rational1 < rational2);
}

bool operator!=(const Rational& rational1, const Rational& rational2) {
    return (rational1 < rational2 || rational2 < rational1);
}

bool operator==(const Rational& rational1, const Rational& rational2) {
    return !(rational1 < rational2 || rational2 < rational1);
}
