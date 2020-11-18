# Modul pro zobrazovani informaci o postupu
# Autor: Bc. Jan Melicharik
# Datum: 24.10.2020

def progress_bar(current_run: int, whole: int, bar_length: int = 20):
    step = round(whole/100)
    if current_run % step == 0:
        marks = (current_run // step) // (100 // bar_length)
        bar = "[" + "=" * marks + " " * (bar_length-marks) + "]"
        print(f"Postup: {bar} {round(100*current_run/whole)}%", end="\r", flush=True)
        
    if current_run + 1 == whole:
        print(f"Postup: [{'='*bar_length}] 100%")
