import matplotlib.pyplot as plt
from collections import deque
import argparse
import os


def calculation(sequence):  # Считает процент GC-состава
    g = sequence.count("G")
    c = sequence.count("C")
    return (g + c) / len(sequence)


def gc_profile(filename, window_size):  # Строит скользящий профиль GC-состава
    profile = []
    window = deque()
    full_seq = []

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip().upper()
            if line == "" or line.startswith(">"):
                continue

            for base in line:
                full_seq.append(base)
                window.append(base)
                if len(window) > window_size:
                    window.popleft()
                if len(window) == window_size:
                    window_str = "".join(window)
                    gc_val = calculation(window_str) * 100
                    profile.append(gc_val)
    total_gc = calculation("".join(full_seq)) * 100
    return total_gc, profile


def main():  # Точка входа в программу (вызывается через pyproject.toml)
    script_dir = os.path.dirname(
        os.path.abspath(__file__)
    )  # Определяем папку, где лежит этот файл (src/gc_profile/)
    my_fasta = os.path.join(
        script_dir, "..", "..", "examples", "sequence.fasta"
    )  # Путь к примеру файла: на 2 уровня выше лежит папка examples/

    parser = argparse.ArgumentParser(description="Построение GC-профиля ДНК")
    parser.add_argument("-i", "--input", default=my_fasta, help="путь к fasta файлу")
    parser.add_argument("-w", "--window", type=int, default=100, help="размер окна")
    parser.add_argument(
        "-o", "--output", default=None, help="имя файла для сохранения графика"
    )
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print("Ошибка: файл не найден.")
        exit()

    print("Файл:", args.input)
    print("Размер окна:", args.window)

    gc_total, profile = gc_profile(args.input, args.window)

    print(f"Общий GC-состав: {gc_total:.2f}%")
    print(f"Точек в профиле: {len(profile)}")

    plt.figure(figsize=(12, 5))
    plt.plot(profile, color="blue", label="GC-профиль по окнам")
    plt.axhline(
        y=gc_total, color="red", linestyle="--", label=f"Средний GC: {gc_total:.2f}%"
    )
    plt.title("GC-профиль последовательности ДНК\n(курсовая работа, Носикова М.С.)")
    plt.xlabel("Позиция (номер окна)")
    plt.ylabel("GC-состав, %")
    plt.ylim(0, 100)
    plt.legend()
    plt.grid(True, alpha=0.3)

    if args.output:
        plt.savefig(args.output)
        print("График сохранен в:", args.output)
    else:
        plt.show()


# Запуск только при прямом вызове файла
if __name__ == "__main__":
    main()
