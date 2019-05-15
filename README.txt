Где скачать PyMOL?

https://pymol.org/2/


Как подключить PyMOL?

В терминале pymol -R
from xmlrpc.client import ServerProxy
pymol = ServerProxy(uri="http://localhost:9123/RPC2") #номер хоста из терминала

Несколько команд:
pymol.fetch("1dn2") #загрузить из PDB белок под этим названием
pymol.delete("1dn2")
pymol.bg_colour("black") #установить цвет фона черный
pymol.show_as("cartoon")
pymol.do('fetch 1bl7')
