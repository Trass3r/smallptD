@del .deps
@del .objs /Q

xfbuild +v smallpt.d +obin\smallptD +xcore +xstd -I..\cl4d -debug -g -unittest && cv2pdb -D2 bin\smallptD.exe