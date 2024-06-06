# Структура проекта

Папка gnuplots — *.plt-файлы — скрипты для gnuplot и *.d-файлы с данными.

Папка plots — готовые графики.

Папка source — исходники 

# Демонстрация работы:

**Вариант 1)** В папке gnuplots собраны скрипты для gnuplot и уже созданные файлы с данными

- gnuplotRossler.plt - Выводит аттрактор Ресслера

- gnuplotRossler_a.plt - Выводит бифуркационную диаграмму для параметра a

- gnuplotRossler_b.plt - Выводит бифуркационную диаграмму для параметра b

- gnuplotRossler_c.plt - Выводит бифуркационную диаграмму для параметра c

**Вариант 2)** В папке source хранятся исходники кода на C++. Они разбиты по папкам:

- Rossler - main.cpp - создает файл Rossler.d с данными для построения через gnuplotRossler.plt

- Rossler_a - main.cpp - создает файл Rossler_a.d с данными для построения через gnuplotRossler_a.plt

- Rossler_b - main.cpp - создает файл Rossler_b.d с данными для построения через gnuplotRossler_b.plt

- Rossler_c - main.cpp - создает файл Rossler_c.d с данными для построения через gnuplotRossler_c.plt

Компиляцию каждого файла main проводить с помощью команды:

`g++ main.c`
