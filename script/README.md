## Скрипты

* [generation.py](generation.py) - для генерации функции распределения тепла в начальный момент на основе [setting.ini](../initial/setting.ini)

* [plotX.py](plotX.py) - для построения сечения вдоль плоскости X графиков распределения в начальный и конечный *(результат работы программы)* момент времени. 

Параметры входящей строки:
1. Название раздела где находятся результаты относительно папки [result](../result)
2. Название файла с результатами.

*пример:*
```
    plotX.py Kirill result.txt
```