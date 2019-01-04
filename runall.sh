#!/bin/bash
FRESCO=frescox #Calling FRESCO
GRACE=xmgrace
destino=fort.4
XSECS2=xsec_states.out # Data from FRESCO's output fort.13 [Xsec to all excited states one by one].
XSECS3=Cross_sections.out # Cross section output for each energy.
fort=fort.13
fort39=fort.39
lista=lista.txt
grace_file=graphs.gr
outfile=outfile.out
aux=aux.dat
XSECS=xsec_completo.out
rm -f $destino; rm -f $XSECS; rm -f $XSECS2; rm -f $outfile; rm -f $aux; 
echo 'Energy Elastic    Absorption   Reaction    Total    (MeV/mb)' >> $aux
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
  echo 'Elastic:' $xel mb
  echo 'Reaction:' $xreac mb
  echo 'Absorption:' $xabs mb
  echo ""
  echo $elab $xel $xabs $xreac $xtot  >> $XSECS
  echo $elab $x_gs_0 >> $XSECS2
  cat $output >> $outfile
  rm -f $output
done < $lista
if [ -d Results ]; then
    echo 'Results of the run moved to "Results" folder'
else   
    mkdir  Results
    echo 'Results of the run moved to "Results" folder'
fi
if [ -e $grace_file ]; then
  $GRACE -batch $grace_file -nosafe -hardcopy 
  mv total.eps Results; mv elastic.eps Results; mv absorption.eps Results; mv $outfile Results
fi
cat $aux $XSECS > $XSECS3
mv $XSECS3 Results; mv $XSECS2 Results
rm -f $destino; rm -f $aux; rm -f $XSECS; rm -f $grace_file
