                   ____                    _   _                         
                  |  _ \   /\        /\   | \ | |                        
                  | |_) | /  \      /  \  |  \| |                        
                  |  _ < / /\ \    / /\ \ | . ` |                        
                  | |_) / ____ \  / ____ \| |\  |                        
  _____       _   |____/_/    \_\/_/    \_\_| \_|         _       _      
 |  __ \     | | | |                                     | |     | |     
 | |__) |   _| |_| |__   ___  _ __    _ __ ___   ___   __| |_   _| | ___ 
 |  ___/ | | | __| '_ \ / _ \| '_ \  | '_ ` _ \ / _ \ / _` | | | | |/ _ \
 | |   | |_| | |_| | | | (_) | | | | | | | | | | (_) | (_| | |_| | |  __/
 |_|    \__, |\__|_| |_|\___/|_| |_| |_| |_| |_|\___/ \__,_|\__,_|_|\___|
         __/ |                                                           
        |___/                                                            

Projekt: Nástroje a techniky bayesiánské ekonometrie a statistiky (MUNI/FR/1108/2019)
Autor: Bc. Jan Melichařík (451747)

Sada skriptů a pomocných funkcí v jazyce Python kopírujících postupy skriptů ze cvičení (program MATLAB).

Naposledy aktualizováno: <<DATUM ZDE>>

--------------------------------------------------------------------------

Pro rychlejší zpřístupnění skriptů je součástí projektu i virtuální prostředí pro jazyk Python, které obsahuje veškeré potřebné knihovny pro spuštění jednotlivých skriptů.

Aby byly skripty funkční, je třeba nastavit cestu v terminálu do složky jednotlivých cvičení. Datové soubory jsou v těchto složkách uloženy také a pokud nebudou skripty spuštěny z této složky, může dojít k problémům s načítáním dat.

Poznámka 1: V Pythonu může dojít k problému, pokud se budete snažit vypočítat mocninu negativního čísla na reálné číslo, tak díky přednosti mocniny před "znaménkem" čísla můžete získat hodnotu "nan" (toto se děje při umocňování matic). Tento typ se dá kontrolovat pomocí funkce numpy.isnan().

