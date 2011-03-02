@del .deps
@del .objs /Q

xfbuild +v smallpt.d +obin\smallptD +xcore +xstd -I..\cl4d -release -O -inline
pause