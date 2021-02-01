# Bayesiánská analýza                                      

Projekt: Nástroje a techniky bayesiánské ekonometrie a statistiky (MUNI/FR/1108/2019)

Autor: Bc. Jan Melichařík (451747@mail.muni.cz)

Sada skriptů a pomocných funkcí v jazyce Python (verze 3.9.1) kopírujících postupy skriptů ze cvičení (program MATLAB).

Naposledy aktualizováno: 1.2.2021

## Obecné informace

Pro rychlejší zpřístupnění skriptů je součástí projektu i virtuální prostředí pro jazyk Python, které obsahuje veškeré potřebné knihovny pro spuštění jednotlivých skriptů.

Aby byly skripty funkční, je třeba nastavit cestu v terminálu do složky jednotlivých cvičení. Datové soubory jsou v těchto složkách uloženy také a pokud nebudou skripty spuštěny z této složky, může dojít k problémům s načítáním dat.

### Použité knihovny

Standardní:
- `math`
- `sys`

Rozšířené:
- `matplotlib`
- `numpy`
- `pandas`
- `scipy`
- `statistics`
- `tabulate`

Pomocné funkce převzaté od doc. Daniela Němce a uloženy v knihovně support:
- `gamm_rnd_koop` podle `gamm_rnd_Koop_2.m`
- `geweke`
- `lik_ces`
- `log_post_ces`
- `my_boot`
- `my_nlrm`
- `norm_rnd`
- `prior_ces`
- `progress_bar`
- ...

## Poznámky

Poznámka 1: V Pythonu může dojít k problému, pokud se budete snažit vypočítat mocninu negativního čísla na reálné číslo, tak díky přednosti mocniny před "znaménkem" čísla můžete získat hodnotu "nan" (toto se děje při umocňování matic). Tento typ se dá kontrolovat pomocí funkce numpy.isnan().


