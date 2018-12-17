#!/bin/bash
FRESCO=frescox
destino=fort.4
XSECS=xsec_completo.out
XSECS2=xsec_estados.out
fort=fort.13
fort39=fort.39
rm -f $destino
rm -f $XSECS
rm -f $XSECS2
while read input origen; do
 base=${input%%.*}
 output=$base.out
 ln -fs $origen $destino
#  echo $input $output
  elab=`cat  $input | egrep elab | awk '{print $2 }'`
  echo Running energy: $elab
  
  $FRESCO  < $input  > $output 

 xel=`grep 0.0000 $fort  | awk '{print $9}'`
 xtot=`grep 0.0000 $fort | awk '{print $8}'` 
 xabs=`grep -v NaN $fort39 | awk '{print $2}'`
 xreac=`grep 0.0000 $fort | awk 'FNR ==2 {print $7}'`
 x_gs_0=`grep  0.5  $fort  | awk '{print $7}'`
  echo 'Elastic:' $xel ' Reac:' $xreac  ' Abs:' $xabs #'0+gs: ' $x_gs_0
  echo $elab $xel $xabs $xreac $xtot  >> $XSECS
  echo $elab $x_gs_0 >> $XSECS2
#  mv fort.56
  
done < lista.txt
rm -f $destino
#Uncomment next line to generate C.S graphs using GRACE in batch mode [script graphs.gr is attached]
#xmgrace -batch graphs.gr -nosafe -hardcopy -log xy




