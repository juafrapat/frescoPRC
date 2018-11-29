#!/bin/bash
FRESCO=frescox
destino=fort.4
XSECS=xsec_completo.out
XSECS2=xsec_estados.out
fort=fort.13
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

 xel=`grep ELASTIC $output  | awk '{print $6}'`
 xtot=`grep TOTAL  $output | awk '{print $6}'` 
 xabs=`grep ABSORBTION $output | awk '{print $7}'`
 xreac=`grep REACTION $output | awk '{print $6}'`
 x_gs_0=`grep  0.5  $fort  | awk '{print $7}'`
  echo 'Elastic:' $xel ' Reac:' $xreac  ' Abs:' $xabs #'0+gs: ' $x_gs_0
  echo $elab $xel $xabs $xreac $xtot  >> $XSECS
  echo $elab $x_gs_0 >> $XSECS2
#  mv fort.56
  
done < lista.txt






