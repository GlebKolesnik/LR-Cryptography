import math
import random


def phi(n):
    """Функция Эйлера."""
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result


def mu(n):
    """Функция Мебиуса."""
    if n == 1:
        return 1
    prime_factors = 0
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            prime_factors += 1
            if n % p == 0:
                return 0
        p += 1
    if n > 1:
        prime_factors += 1
    return -1 if prime_factors % 2 == 1 else 1


def lcm(a, b):
    """Наименьшее общее кратное."""
    return abs(a * b) // math.gcd(a, b)


def lcm_list(numbers):
    """Наименьшее общее кратное для списка чисел."""
    lcm_result = 1
    for num in numbers:
        lcm_result = lcm(lcm_result, num)
    return lcm_result


# Тесты
test_number = 12
test_numbers_list = [12, 15, 20]

phi(test_number), mu(test_number), lcm_list(test_numbers_list)


def extended_gcd(a, b):
    """Расширенный алгоритм Евклида."""
    if b == 0:
        return (a, 1, 0)
    else:
        d, x, y = extended_gcd(b, a % b)
        return (d, y, x - (a // b) * y)


def chinese_remainder_theorem(n, a):
    """Китайская теорема остатков."""
    prod = 1
    for ni in n:
        prod *= ni

    result = 0
    for ni, ai in zip(n, a):
        p = prod // ni
        _, x, _ = extended_gcd(p, ni)
        result += ai * x * p

    return result % prod


# Тесты
test_moduli = [3, 5, 7]
test_remainders = [2, 3, 2]

chinese_remainder_theorem(test_moduli, test_remainders)


def legendre_symbol(a, p):
    """Символ Лежандра."""
    if a % p == 0:
        return 0
    elif pow(a, (p - 1) // 2, p) == 1:
        return 1
    else:
        return -1


def jacobi_symbol(a, n):
    """Символ Якоби."""
    if n % 2 == 0 or n <= 0:
        raise ValueError("n must be odd and positive.")
    a %= n
    if a == 0:
        return 0
    if a == 1:
        return 1
    if a == 2:
        if n % 8 in [3, 5]:
            return -1
        else:
            return 1
    if a == n - 1:
        if n % 4 == 1:
            return 1
        else:
            return -1
    if math.gcd(a, n) > 1:
        return 0
    return jacobi_symbol(n % a, a) * (-1 if (a % 4 == 3 and n % 4 == 3) else 1)


# Тесты
test_a = 5
test_p = 7
test_n = 15

legendre_result = legendre_symbol(test_a, test_p)
jacobi_result = jacobi_symbol(test_a, test_n)

legendre_result, jacobi_result


def pollard_rho(n, max_iterations=10000):
    """Ро-алгоритм Полларда для факторизации."""
    if n % 2 == 0:
        return 2
    x = 2;
    y = 2;
    d = 1
    f = lambda x: (x ** 2 + 1) % n
    i = 0
    while d == 1 and i < max_iterations:
        x = f(x)
        y = f(f(y))
        d = math.gcd(abs(x - y), n)
        i += 1
    return d


# Тест
test_number = 8051
pollard_rho(test_number)


def baby_step_giant_step(g, h, p):
    """Алгоритм 'большой шаг – малый шаг' для дискретного логарифма."""
    m = int(p ** 0.5) + 1
    lookup = {}

    # Большой шаг
    for j in range(m):
        value = pow(g, j, p)
        lookup[value] = j

    # Находим инверсию g по модулю p
    inv_g = pow(g, -m, p)
    current_h = h

    # Малый шаг
    for i in range(m):
        if current_h in lookup:
            return i * m + lookup[current_h]
        else:
            current_h = current_h * inv_g % p

    raise ValueError("Дискретный логарифм не найден.")


# Тест
test_g = 2
test_h = 22
test_p = 29

baby_step_giant_step(test_g, test_h, test_p)


def jacobi_symbol(a, n):
    """Улучшенный символ Якоби."""
    if n == 1:
        return 1
    if a == 0:
        return 0
    if a == 1:
        return 1
    if a == 2:
        if n % 8 in [3, 5]:
            return -1
        else:
            return 1
    if a % 2 == 0:
        return jacobi_symbol(2, n) * jacobi_symbol(a // 2, n)

    if a % 4 == 3 and n % 4 == 3:
        return -jacobi_symbol(n % a, a)
    else:
        return jacobi_symbol(n % a, a)


def tonelli_shanks(a, p):
    """Алгоритм Чипполы нахождения дискретного квадратного корня."""

    if legendre_symbol(a, p) != 1:
        return None  # Квадратного корня не существует

    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1

    if s == 1:
        return pow(a, (p + 1) // 4, p)

    # Находим з какое не является квадратичным вычетом по модулю p
    z = 1
    while legendre_symbol(z, p) != -1:
        z += 1

    c = pow(z, q, p)
    r = pow(a, (q + 1) // 2, p)
    t = pow(a, q, p)
    m = s
    t_2 = 0

    while (t - 1) % p != 0:
        t_2 = (t * t) % p
        for i in range(1, m):
            if (t_2 - 1) % p == 0:
                break
            t_2 = (t_2 * t_2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i

    return r


# Тест
test_a = 10
test_p = 29

tonelli_shanks(test_a, test_p)


def solovay_strassen_test(n, k=5):
    """Алгоритм Соловея-Штрассена."""
    if n == 2 or n == 3:
        return True
    if n < 2 or n % 2 == 0:
        return False

    for _ in range(k):
        a = random.randint(2, n - 1)
        x = jacobi_symbol(a, n)
        y = pow(a, (n - 1) // 2, n)
        if x == 0 or y != x % n:
            return False
    return True


# Тест
test_number = 29
print(solovay_strassen_test(test_number))


# Функция для расширенного алгоритма Евклида
def extended_gcd(a, b):
    if a == 0:
        return b, 0, 1
    gcd, x1, y1 = extended_gcd(b % a, a)
    x = y1 - (b // a) * x1
    y = x1
    return gcd, x, y


# Функция для теста Соловея-Штрассена
def solovay_strassen_test(n, k=5):
    def jacobi_symbol(a, n):
        if n == 1:
            return 1
        if a == 0:
            return 0
        if a == 1:
            return 1
        if a == 2:
            if n % 8 in [3, 5]:
                return -1
            else:
                return 1
        if a % 2 == 0:
            return jacobi_symbol(2, n) * jacobi_symbol(a // 2, n)
        if a % 4 == 3 and n % 4 == 3:
            return -jacobi_symbol(n % a, a)
        else:
            return jacobi_symbol(n % a, a)

    if n == 2 or n == 3:
        return True
    if n < 2 or n % 2 == 0:
        return False

    for _ in range(k):
        a = random.randint(2, n - 1)
        x = jacobi_symbol(a, n)
        y = pow(a, (n - 1) // 2, n)
        if x == 0 or y != x % n:
            return False
    return True


# Функция для генерации случайного простого числа
def random_prime(bits):
    while True:
        number = random.getrandbits(bits)
        if solovay_strassen_test(number):
            return number


# Функция для генерации ключей RSA
def generate_rsa_keys(bits=256):
    p = random_prime(bits // 2)
    q = random_prime(bits // 2)
    n = p * q
    phi_n = (p - 1) * (q - 1)

    e = random.randint(2, phi_n - 1)
    while math.gcd(e, phi_n) != 1:
        e = random.randint(2, phi_n - 1)

    _, d, _ = extended_gcd(e, phi_n)
    d %= phi_n
    if d < 0:
        d += phi_n

    return (e, n), (d, n)


# Функция для шифрования сообщения с помощью RSA
def rsa_encrypt(message, public_key):
    e, n = public_key
    return pow(message, e, n)


# Функция для дешифрования шифртекста с помощью RSA
def rsa_decrypt(ciphertext, private_key):
    d, n = private_key
    return pow(ciphertext, d, n)


# Криптосистема Эль-Гамаля
def find_generator(p):
    """Находит генератор группы Z_p^*."""
    factors = [2, (p - 1) // 2]
    for g in range(2, p):
        if all(pow(g, f, p) != 1 for f in factors):
            return g
    raise ValueError("Генератор не найден.")


def generate_elgamal_keys(bits=256):
    """Генерация ключей для Эль-Гамаля."""
    p = random_prime(bits)
    g = find_generator(p)
    x = random.randint(1, p - 2)
    h = pow(g, x, p)

    return (p, g, h), x


def elgamal_encrypt(message, public_key):
    """Шифрование сообщения с помощью Эль-Гамаля."""
    p, g, h = public_key
    y = random.randint(1, p - 2)
    c1 = pow(g, y, p)
    c2 = (message * pow(h, y, p)) % p
    return (c1, c2)


def elgamal_decrypt(ciphertext, private_key, public_key):
    """Дешифрование шифртекста с помощью Эль-Гамаля."""
    p, _, _ = public_key
    x = private_key
    c1, c2 = ciphertext
    s = pow(c1, x, p)
    m = (c2 * pow(s, -1, p)) % p
    return m


# Генерация ключей и тестирование Эль-Гамаля
public_key, private_key = generate_elgamal_keys(bits=256)
test_message_elgamal = 123456789
ciphertext_elgamal = elgamal_encrypt(test_message_elgamal, public_key)
decrypted_message_elgamal = elgamal_decrypt(ciphertext_elgamal, private_key, public_key)

test_message_elgamal, decrypted_message_elgamal


class EllipticCurve:
    """Класс для работы с эллиптической кривой y^2 = x^3 + ax + b (mod p)."""

    def __init__(self, a, b, p):
        self.a = a
        self.b = b
        self.p = p

    def is_valid_point(self, x, y):
        """Проверяет, является ли точка (x, y) допустимой на кривой."""
        return (y * y - x * x * x - self.a * x - self.b) % self.p == 0

    def add(self, P, Q):
        """Складывает две точки P и Q на кривой."""
        if P is None:  # Нейтральный элемент
            return Q
        if Q is None:  # Нейтральный элемент
            return P
        x1, y1 = P
        x2, y2 = Q
        if P == Q:  # Дублирование точки
            if y1 == 0:
                return None
            l = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, self.p)
        else:
            if x1 == x2:
                return None
            l = (y2 - y1) * pow(x2 - x1, -1, self.p)
        x3 = (l * l - x1 - x2) % self.p
        y3 = (l * (x1 - x3) - y1) % self.p
        return x3, y3

    def multiply(self, P, n):
        """Умножает точку P на целое число n."""
        R = None
        m2 = P
        while n:
            if n & 1:
                R = self.add(R, m2)
            m2 = self.add(m2, m2)
            n >>= 1
        return R


def generate_elgamal_ec_keys(curve, G, n):
    """Генерация ключей для Эль-Гамаля на основе эллиптических кривых."""
    x = random.randint(1, n - 1)
    H = curve.multiply(G, x)

    return (curve, G, H, n), x


def elgamal_ec_encrypt(message_point, public_key):
    """Шифрование сообщения (точки) с помощью Эль-Гамаля на основе эллиптических кривых."""
    curve, G, H, n = public_key
    y = random.randint(1, n - 1)
    C1 = curve.multiply(G, y)
    S = curve.multiply(H, y)
    C2 = curve.add(message_point, S)

    return (C1, C2)


def elgamal_ec_decrypt(ciphertext, private_key, public_key):
    """Дешифрование шифртекста с помощью Эль-Гамаля на основе эллиптических кривых."""
    curve, _, _, _ = public_key
    x = private_key
    C1, C2 = ciphertext
    S = curve.multiply(C1, x)
    S_inv = (S[0], -S[1] % curve.p)
    M = curve.add(C2, S_inv)

    return M


def main():
    # 1) Вычисление функций Эйлера и Мебиуса
    n = 30
    print(f"Функция Эйлера для {n}: {phi(n)}")
    print(f"Функция Мебиуса для {n}: {mu(n)}")

    # 2) Китайская теорема об остатках
    mods = [3, 5, 7]
    residues = [2, 3, 2]
    print(f"Решение системы сравнений {residues} по модулям {mods}: {chinese_remainder_theorem(residues, mods)}")

    # 3) Вычисление символов Лежандра и Якоби
    a, p = 5, 11
    print(f"Символ Лежандра ({a}/{p}) = {legendre_symbol(a, p)}")

    # 4) ро-алгоритм Полларда для факторизации
    n = 8051
    print(f"Один из делителей {n} (ро-алгоритм Полларда): {pollard_rho(n)}")

    # 5) ро-алгоритм Полларда для дискретного логарифма
    base, target, modulo = 2, 22, 29
    print(
        f"Дискретный логарифм {base}^x = {target} (mod {modulo}): {baby_step_giant_step(base, target, modulo)}")

    # 6) Алгоритм чипполов
    n, p = 10, 13
    print(f"Дискретный квадратный корень {n} в GF({p}): {tonelli_shanks(n, p)}")

    # 7) Алгоритм Соловея-Штрассена
    n = 561
    print(f"Число {n} простое (по алгоритму Соловея-Штрассена)? {'Да' if solovay_strassen_test(n) else 'Нет'}")

    # 8) RSA
    public_key, private_key = generate_rsa_keys(bits=256)
    test_message_rsa = 123456789
    ciphertext_rsa = rsa_encrypt(test_message_rsa, public_key)
    decrypted_message_rsa = rsa_decrypt(ciphertext_rsa, private_key)
    print(f"RSA: Отправленное сообщение: {test_message_rsa}, Расшифрованное сообщение: {decrypted_message_rsa}")

    # 9) Эль-Гамаль на эллиптических кривых
    curve = EllipticCurve(a=2, b=2, p=17)
    G = (5, 1)
    n = 19
    public_key_ec, private_key_ec = generate_elgamal_ec_keys(curve, G, n)
    test_message_point = (9, 16)
    ciphertext_ec = elgamal_ec_encrypt(test_message_point, public_key_ec)
    decrypted_message_point = elgamal_ec_decrypt(ciphertext_ec, private_key_ec, public_key_ec)
    # Вывод результата функции нахождения НОК
    numbers = [12, 15, 20]
    print(f"Наименьшее общее кратное для чисел {numbers}: {lcm_list(numbers)}")

    # Вывод результата функции вычисления символа Якоби
    a, n = 5, 11
    print(f"Символ Якоби ({a}/{n}) = {jacobi_symbol(a, n)}")
    print(
        f"Эль-Гамаль на эллиптических кривых: Отправленное сообщение: {test_message_point}, Расшифрованное сообщение: {decrypted_message_point}")


if __name__ == "__main__":
    main()
