#!/bin/bash
FRESCO=frescox #Calling FRESCO
destino=fort.4
XSECS=xsec_completo.out # Cross section output for each energy.
XSECS2=xsec_estados.out # Data from FRESCO's output fort.13 [Xsec to all excited states one by one].
fort=fort.13
fort39=fort.39
lista=lista.txt
rm -f $destino
rm -f $XSECS
rm -f $XSECS2
echo 'Energy Elastic    Absorption   Reaction    Total    (MeV/mb)' >> $XSECS
while read input origen; do
 base=${input%%.*}
 output=$base.out
# Comment next line if you will make a run without any excited band (only G.S band).
 ln -fs $origen $destino
  elab=`cat  $input | egrep elab | awk '{print $2 }'`
  echo Running energy: $elab MeV
  
  $FRESCO  < $input  > $output 

 xel=`grep 0.0000 $fort  | awk '{print $9}'`
 xtot=`grep 0.0000 $fort | awk '{print $8}'` 
 xabs=`grep -v NaN $fort39 | awk '{print $2}'`
 xreac=`grep 0.0000 $fort | awk 'FNR ==2 {print $7}'`
 x_gs_0=`grep  0.5  $fort  | awk '{print $7}'`
  echo ""
  echo "Elastic:" $xel mb 
  echo 'Reaction:' $xreac mb
  echo 'Absorption:' $xabs mb
  echo ""
  echo $elab $xel $xabs $xreac $xtot  >> $XSECS
  echo $elab $x_gs_0 >> $XSECS2

  
done < $lista 
rm -f $destino
#Uncomment next line to generate C.S graphs using GRACE in batch mode [script graphs.gr is attached]
#xmgrace -batch graphs.gr -nosafe -hardcopy -log xy




