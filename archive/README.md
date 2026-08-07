# Archivo — archivos en desuso

Archivos que ya no forman parte del pipeline activo (`01_copulas.R`), movidos
acá para no perder historial pero sacarlos del flujo principal del repo.

| Archivo | Path original | Por qué se archivó |
|---|---|---|
| `01_script_maestro.R` | `01_script_maestro.R` | Prototipo exploratorio de una sola estación: hardcodea `estacion.usar` y carga un `.rda` desde una ruta local de otro colaborador (`/home/dbonhaure/RStudioProjects/Copulas/...`). Ningún otro script lo referencia. Superado funcionalmente por `01_copulas.R`, que es el pipeline de producción dirigido por YAML. |
| `funciones_ajuste.xls` | `data/input/funciones_ajuste.xls` | Config de ajuste univariado leída vía `readxl::read_excel(...)`. Esa lectura está comentada en el código (solo la referenciaba, comentada, `01_script_maestro.R`) — fue reemplazada por la tabla `configuracion.ajuste.univariado` embebida en `parametros_copulas.yml`. |
| `funciones_ajuste_copulas.xls` | `data/input/funciones_ajuste_copulas.xls` | Ídem anterior, para la configuración de ajuste de cópulas (`configuracion.ajuste.copula` en `parametros_copulas.yml`). |
