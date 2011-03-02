@del .deps
@del .objs /Q

xfbuild +v smallpt.d +obin\smallptD +xcore +xstd -I..\cl4d -debug -g -L/STACK:20971520 && cv2pdb -D2 bin\smallptD.exe
pause