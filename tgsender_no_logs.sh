#!/bin/bash
chat_Id="639261746"
#chat_IdRinat="422366742"
#chat_IdSegAbr="653206413"
#chat_IDTea="1001581290556"
#chat_IdBelok="-1001337146814"
botToken="1992203014:AAGXCU5ta31M-R10axejbBtxRJd0L1PNOow"
curdir="/media/ElissarDisk/ADASTRA/logs/"
texti='Оно прогналось, пожалуйста, загляните в логи'
first='parallel_log'
second='whole_log'
#botFather
#curl -s -X POST https://api.telegram.org/bot$botToken/sendMessage -d chat_id=$chatId -d text="$texti"

#echo sending $curdir/$1
#curl -F chat_id=$chatId -F document=@$curdir/$1 https://api.telegram.org/bot$botToken/sendDocument

#echo sending $curdir/$2
#curl -F chat_id=$chatId -F document=@$curdir/$2 https://api.telegram.org/bot$botToken/sendDocument

#########################
curl -s -X POST https://api.telegram.org/bot$botToken/sendMessage -d chat_id=$chat_Id -d text="$texti"
curl -s -X POST https://api.telegram.org/bot$botToken/sendMessage -d chat_id=$chat_IdBelok -d text="$texti"
echo sending $curdir/$first
curl -F chat_id=$chat_IdBelok -F document=$curdir/$first https://api.telegram.org/bot$botToken/sendDocument
echo sending $curdir/$second
curl -F chat_id=$chat_IdBelok -F document=@$curdir/$second https://api.telegram.org/bot$botToken/sendDocument