### Функция Эйлера и Мебиуса (`phi` и `mu`)

- **Функция Эйлера (Фи-функция)**: Данная функция вычисляет количество чисел, меньших \( n \) и взаимно простых с \( n \).
- **Функция Мебиуса**: Возвращает `-1`, если \( n \) имеет нечетное число простых делителей (и не имеет кратных простых делителей), `1` в противном случае, и `0`, если \( n \) имеет кратный простой делитель.

### Наименьшее общее кратное (`lcm` и `lcm_list`)

- **Наименьшее общее кратное (НОК)**: Функция вычисляет НОК для списка чисел, используя стандартный алгоритм, основанный на НОД (наибольший общий делитель).

### Китайская теорема об остатках (`chinese_remainder_theorem`)

- Данная функция решает систему линейных сравнений и возвращает общее решение. Она использует расширенный алгоритм Евклида для нахождения обратных элементов.

### Символы Лежандра и Якоби (`legendre_symbol` и `jacobi_symbol`)

- **Символ Лежандра**: Для простого числа \( p \) и целого числа \( a \), символ Лежандра определяет, является ли \( a \) квадратичным вычетом по модулю \( p \).
- **Символ Якоби**: Обобщение символа Лежандра для любого нечетного натурального числа \( n \).

### Ро-алгоритм Полларда для факторизации (`pollard_rho`)

- Этот алгоритм позволяет находить нетривиальные делители составного числа. Он может быть особенно эффективным для чисел с небольшими простыми делителями.

### Алгоритм нахождения дискретного логарифма (`baby_step_giant_step`)

- Данный алгоритм использует метод "шаг младенца - шаг великана" для вычисления дискретных логарифмов в конечных полях.

### Алгоритм Чиполли нахождения дискретного квадратного корня (`tonelli_shanks`)

- Этот алгоритм находит квадратный корень числа \( a \) по модулю простого числа \( p \), если таковой существует.

### Алгоритм Соловея-Штрассена (`solovay_strassen_test`)

- Это вероятностный тест для проверки чисел на простоту. Он использует свойства символа Якоби и быстро выявляет составные числа.

### Криптосистема RSA (`generate_rsa_keys`, `rsa_encrypt`, `rsa_decrypt`)

- **RSA** - это криптосистема с открытым ключом, которая использует свойства больших простых чисел и функции Эйлера. В вашем коде реализованы функции для генерации ключей, шифрования и дешифрования сообщений.

### Криптосистема Эль-Гамаля на эллиптических кривых (`generate_elgamal_ec_keys`, `elgamal_ec_encrypt`, `elgamal_ec_decrypt`)

- Это обобщение традиционной криптосистемы Эль-Гамаля на эллиптические кривые. Она предоставляет улучшенную безопасность при меньших размерах ключей.
