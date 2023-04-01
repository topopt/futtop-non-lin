type localMatrix = {xx: f64, xy: f64, xz: f64, yx: f64, yy: f64, yz: f64, zx: f64, zy: f64, zz: f64}

let getke_l0(recv :i32, send :i32)  :localMatrix = 
match (recv,send)
case (0,0) -> {xx=(0.2350427350427349848516201),xy=(-0.0801282051282050933327383),xz=(0.0801282051282050933327383),yx=(-0.0801282051282050933327383),yy=(0.2350427350427349848516201),yz=(-0.0801282051282050933327383),zx=(0.0801282051282050933327383),zy=(-0.0801282051282050933327383),zz=(0.2350427350427349570960445)}
case (0,1) -> {xx=(-0.1068376068376068188658934),xy=(-0.0160256410256410172787689),xz=(0.0160256410256410172787689),yx=(0.0160256410256410138093219),yy=(0.0534188034188033816773711),yz=(-0.0400641025641025397274753),zx=(-0.0160256410256410172787689),zy=(-0.0400641025641025397274753),zz=(0.0534188034188033816773711)}
case (0,2) -> {xx=(-0.0854700854700854301126967),xy=(0.0801282051282050933327383),xz=(0.0080128205128205121088314),yx=(0.0801282051282050794549505),yy=(-0.0854700854700854301126967),yz=(-0.0080128205128205138435549),zx=(-0.0080128205128205034352140),zy=(0.0080128205128205034352140),zz=(-0.0053418803418803411167670)}
case (0,3) -> {xx=(0.0534188034188033816773711),xy=(0.0160256410256410172787689),xz=(0.0400641025641025327885814),yx=(-0.0160256410256410172787689),yy=(-0.1068376068376068049881056),yz=(-0.0160256410256410172787689),zx=(0.0400641025641025258496875),zy=(0.0160256410256410138093219),zz=(0.0534188034188033886162650)}
case (0,4) -> {xx=(0.0534188034188033816773711),xy=(-0.0400641025641025327885814),xz=(-0.0160256410256410172787689),yx=(-0.0400641025641025327885814),yy=(0.0534188034188033816773711),yz=(0.0160256410256410207482158),zx=(0.0160256410256410172787689),zy=(-0.0160256410256410172787689),zz=(-0.1068376068376067911103178)}
case (0,5) -> {xx=(-0.0854700854700854301126967),xy=(-0.0080128205128205086393844),xz=(-0.0801282051282050933327383),yx=(0.0080128205128205086393844),yy=(-0.0053418803418803480556609),yz=(0.0080128205128205103741079),zx=(-0.0801282051282050933327383),zy=(-0.0080128205128205069046610),zz=(-0.0854700854700854439904845)}
case (0,6) -> {xx=(-0.0587606837606837115184355),xy=(0.0400641025641025397274753),xz=(-0.0400641025641025327885814),yx=(0.0400641025641025397274753),yy=(-0.0587606837606837184573294),yz=(0.0400641025641025397274753),zx=(-0.0400641025641025327885814),zy=(0.0400641025641025397274753),zz=(-0.0587606837606837115184355)}
case (0,7) -> {xx=(-0.0053418803418803480556609),xy=(0.0080128205128205086393844),xz=(-0.0080128205128205051699375),yx=(-0.0080128205128205051699375),yy=(-0.0854700854700854439904845),yz=(0.0801282051282050933327383),zx=(0.0080128205128205069046610),zy=(0.0801282051282050794549505),zz=(-0.0854700854700854301126967)}
case (1,0) -> {xx=(-0.1068376068376068049881056),xy=(0.0160256410256410138093219),xz=(-0.0160256410256410172787689),yx=(-0.0160256410256410172787689),yy=(0.0534188034188033816773711),yz=(-0.0400641025641025397274753),zx=(0.0160256410256410172787689),zy=(-0.0400641025641025327885814),zz=(0.0534188034188033816773711)}
case (1,1) -> {xx=(0.2350427350427349293404689),xy=(0.0801282051282050933327383),xz=(-0.0801282051282050933327383),yx=(0.0801282051282050933327383),yy=(0.2350427350427349293404689),yz=(-0.0801282051282050933327383),zx=(-0.0801282051282050933327383),zy=(-0.0801282051282050933327383),zz=(0.2350427350427349293404689)}
case (1,2) -> {xx=(0.0534188034188033816773711),xy=(-0.0160256410256410172787689),xz=(-0.0400641025641025397274753),yx=(0.0160256410256410172787689),yy=(-0.1068376068376067911103178),yz=(-0.0160256410256410311565567),zx=(-0.0400641025641025327885814),zy=(0.0160256410256410103398750),zz=(0.0534188034188033816773711)}
case (1,3) -> {xx=(-0.0854700854700854301126967),xy=(-0.0801282051282050794549505),xz=(-0.0080128205128205121088314),yx=(-0.0801282051282050933327383),yy=(-0.0854700854700854301126967),yz=(-0.0080128205128205103741079),zx=(0.0080128205128205017004905),zy=(0.0080128205128205034352140),zz=(-0.0053418803418803428514905)}
case (1,4) -> {xx=(-0.0854700854700854301126967),xy=(0.0080128205128205086393844),xz=(0.0801282051282050933327383),yx=(-0.0080128205128205086393844),yy=(-0.0053418803418803454535757),yz=(0.0080128205128205086393844),zx=(0.0801282051282050794549505),zy=(-0.0080128205128205086393844),zz=(-0.0854700854700854301126967)}
case (1,5) -> {xx=(0.0534188034188033747384772),xy=(0.0400641025641025397274753),xz=(0.0160256410256410172787689),yx=(0.0400641025641025397274753),yy=(0.0534188034188033747384772),yz=(0.0160256410256410172787689),zx=(-0.0160256410256410172787689),zy=(-0.0160256410256410172787689),zz=(-0.1068376068376068049881056)}
case (1,6) -> {xx=(-0.0053418803418803506577461),xy=(-0.0080128205128205086393844),xz=(0.0080128205128205086393844),yx=(0.0080128205128205086393844),yy=(-0.0854700854700854301126967),yz=(0.0801282051282050933327383),zx=(-0.0080128205128205069046610),zy=(0.0801282051282050794549505),zz=(-0.0854700854700854301126967)}
case (1,7) -> {xx=(-0.0587606837606837115184355),xy=(-0.0400641025641025327885814),xz=(0.0400641025641025397274753),yx=(-0.0400641025641025327885814),yy=(-0.0587606837606837184573294),yz=(0.0400641025641025397274753),zx=(0.0400641025641025327885814),zy=(0.0400641025641025397274753),zz=(-0.0587606837606837184573294)}
case (2,0) -> {xx=(-0.0854700854700854301126967),xy=(0.0801282051282050794549505),xz=(-0.0080128205128205034352140),yx=(0.0801282051282050933327383),yy=(-0.0854700854700854301126967),yz=(0.0080128205128205051699375),zx=(0.0080128205128205121088314),zy=(-0.0080128205128205103741079),zz=(-0.0053418803418803411167670)}
case (2,1) -> {xx=(0.0534188034188033816773711),xy=(0.0160256410256410207482158),xz=(-0.0400641025641025327885814),yx=(-0.0160256410256410172787689),yy=(-0.1068376068376067911103178),yz=(0.0160256410256410103398750),zx=(-0.0400641025641025397274753),zy=(-0.0160256410256410242176628),zz=(0.0534188034188033816773711)}
case (2,2) -> {xx=(0.2350427350427349570960445),xy=(-0.0801282051282050933327383),xz=(-0.0801282051282051210883139),yx=(-0.0801282051282050933327383),yy=(0.2350427350427349570960445),yz=(0.0801282051282051210883139),zx=(-0.0801282051282051210883139),zy=(0.0801282051282051210883139),zz=(0.2350427350427350403627713)}
case (2,3) -> {xx=(-0.1068376068376068049881056),xy=(-0.0160256410256410172787689),xz=(-0.0160256410256410103398750),yx=(0.0160256410256410172787689),yy=(0.0534188034188033816773711),yz=(0.0400641025641025397274753),zx=(0.0160256410256410276871097),zy=(0.0400641025641025397274753),zz=(0.0534188034188033955551589)}
case (2,4) -> {xx=(-0.0587606837606837184573294),xy=(0.0400641025641025397274753),xz=(0.0400641025641025397274753),yx=(0.0400641025641025397274753),yy=(-0.0587606837606837184573294),yz=(-0.0400641025641025397274753),zx=(0.0400641025641025397274753),zy=(-0.0400641025641025397274753),zz=(-0.0587606837606837253962233)}
case (2,5) -> {xx=(-0.0053418803418803558619166),xy=(0.0080128205128205069046610),xz=(0.0080128205128205086393844),yx=(-0.0080128205128205086393844),yy=(-0.0854700854700854439904845),yz=(-0.0801282051282050933327383),zx=(-0.0080128205128205086393844),zy=(-0.0801282051282051072105261),zz=(-0.0854700854700854439904845)}
case (2,6) -> {xx=(0.0534188034188033677995833),xy=(-0.0400641025641025327885814),xz=(0.0160256410256410138093219),yx=(-0.0400641025641025327885814),yy=(0.0534188034188033816773711),yz=(-0.0160256410256410172787689),zx=(-0.0160256410256410207482158),zy=(0.0160256410256410242176628),zz=(-0.1068376068376068466214690)}
case (2,7) -> {xx=(-0.0854700854700854162349088),xy=(-0.0080128205128205051699375),xz=(0.0801282051282050794549505),yx=(0.0080128205128205086393844),yy=(-0.0053418803418803463209374),yz=(-0.0080128205128205017004905),zx=(0.0801282051282051072105261),zy=(0.0080128205128205121088314),zz=(-0.0854700854700854439904845)}
case (3,0) -> {xx=(0.0534188034188033816773711),xy=(-0.0160256410256410172787689),xz=(0.0400641025641025327885814),yx=(0.0160256410256410172787689),yy=(-0.1068376068376068188658934),yz=(0.0160256410256410172787689),zx=(0.0400641025641025327885814),zy=(-0.0160256410256410172787689),zz=(0.0534188034188033886162650)}
case (3,1) -> {xx=(-0.0854700854700854301126967),xy=(-0.0801282051282050933327383),xz=(0.0080128205128205017004905),yx=(-0.0801282051282050794549505),yy=(-0.0854700854700854301126967),yz=(0.0080128205128205034352140),zx=(-0.0080128205128205103741079),zy=(-0.0080128205128205103741079),zz=(-0.0053418803418803437188522)}
case (3,2) -> {xx=(-0.1068376068376068049881056),xy=(0.0160256410256410172787689),xz=(0.0160256410256410276871097),yx=(-0.0160256410256410172787689),yy=(0.0534188034188033816773711),yz=(0.0400641025641025397274753),zx=(-0.0160256410256410068704280),zy=(0.0400641025641025397274753),zz=(0.0534188034188034094329467)}
case (3,3) -> {xx=(0.2350427350427349848516201),xy=(0.0801282051282050933327383),xz=(0.0801282051282050933327383),yx=(0.0801282051282050933327383),yy=(0.2350427350427349848516201),yz=(0.0801282051282051072105261),zx=(0.0801282051282050933327383),zy=(0.0801282051282051072105261),zz=(0.2350427350427349848516201)}
case (3,4) -> {xx=(-0.0053418803418803575966400),xy=(-0.0080128205128205051699375),xz=(-0.0080128205128205069046610),yx=(0.0080128205128205086393844),yy=(-0.0854700854700854439904845),yz=(-0.0801282051282051072105261),zx=(0.0080128205128205086393844),zy=(-0.0801282051282050933327383),zz=(-0.0854700854700854439904845)}
case (3,5) -> {xx=(-0.0587606837606837115184355),xy=(-0.0400641025641025327885814),xz=(-0.0400641025641025397274753),yx=(-0.0400641025641025327885814),yy=(-0.0587606837606837184573294),yz=(-0.0400641025641025397274753),zx=(-0.0400641025641025397274753),zy=(-0.0400641025641025397274753),zz=(-0.0587606837606837184573294)}
case (3,6) -> {xx=(-0.0854700854700854162349088),xy=(0.0080128205128205051699375),xz=(-0.0801282051282050933327383),yx=(-0.0080128205128205051699375),yy=(-0.0053418803418803437188522),yz=(-0.0080128205128205051699375),zx=(-0.0801282051282050794549505),zy=(0.0080128205128205121088314),zz=(-0.0854700854700854439904845)}
case (3,7) -> {xx=(0.0534188034188033608606894),xy=(0.0400641025641025327885814),xz=(-0.0160256410256410103398750),yx=(0.0400641025641025327885814),yy=(0.0534188034188033816773711),yz=(-0.0160256410256410103398750),zx=(0.0160256410256410207482158),zy=(0.0160256410256410242176628),zz=(-0.1068376068376068188658934)}
case (4,0) -> {xx=(0.0534188034188033816773711),xy=(-0.0400641025641025327885814),xz=(0.0160256410256410172787689),yx=(-0.0400641025641025327885814),yy=(0.0534188034188033816773711),yz=(-0.0160256410256410172787689),zx=(-0.0160256410256410172787689),zy=(0.0160256410256410207482158),zz=(-0.1068376068376067911103178)}
case (4,1) -> {xx=(-0.0854700854700854301126967),xy=(-0.0080128205128205086393844),xz=(0.0801282051282050933327383),yx=(0.0080128205128205086393844),yy=(-0.0053418803418803454535757),yz=(-0.0080128205128205086393844),zx=(0.0801282051282050933327383),zy=(0.0080128205128205069046610),zz=(-0.0854700854700854301126967)}
case (4,2) -> {xx=(-0.0587606837606837184573294),xy=(0.0400641025641025397274753),xz=(0.0400641025641025397274753),yx=(0.0400641025641025397274753),yy=(-0.0587606837606837184573294),yz=(-0.0400641025641025397274753),zx=(0.0400641025641025397274753),zy=(-0.0400641025641025397274753),zz=(-0.0587606837606837253962233)}
case (4,3) -> {xx=(-0.0053418803418803575966400),xy=(0.0080128205128205069046610),xz=(0.0080128205128205069046610),yx=(-0.0080128205128205086393844),yy=(-0.0854700854700854439904845),yz=(-0.0801282051282050933327383),zx=(-0.0080128205128205069046610),zy=(-0.0801282051282051072105261),zz=(-0.0854700854700854439904845)}
case (4,4) -> {xx=(0.2350427350427349848516201),xy=(-0.0801282051282050933327383),xz=(-0.0801282051282050933327383),yx=(-0.0801282051282050933327383),yy=(0.2350427350427349848516201),yz=(0.0801282051282050933327383),zx=(-0.0801282051282050933327383),zy=(0.0801282051282050933327383),zz=(0.2350427350427349848516201)}
case (4,5) -> {xx=(-0.1068376068376067911103178),xy=(-0.0160256410256410138093219),xz=(-0.0160256410256410207482158),yx=(0.0160256410256410172787689),yy=(0.0534188034188033955551589),yz=(0.0400641025641025466663692),zx=(0.0160256410256410103398750),zy=(0.0400641025641025466663692),zz=(0.0534188034188033886162650)}
case (4,6) -> {xx=(-0.0854700854700854301126967),xy=(0.0801282051282050933327383),xz=(-0.0080128205128205034352140),yx=(0.0801282051282050933327383),yy=(-0.0854700854700854439904845),yz=(0.0080128205128205086393844),zx=(0.0080128205128205103741079),zy=(-0.0080128205128205086393844),zz=(-0.0053418803418803428514905)}
case (4,7) -> {xx=(0.0534188034188033677995833),xy=(0.0160256410256410172787689),xz=(-0.0400641025641025327885814),yx=(-0.0160256410256410172787689),yy=(-0.1068376068376068188658934),yz=(0.0160256410256410172787689),zx=(-0.0400641025641025397274753),zy=(-0.0160256410256410172787689),zz=(0.0534188034188033955551589)}
case (5,0) -> {xx=(-0.0854700854700854301126967),xy=(0.0080128205128205069046610),xz=(-0.0801282051282050933327383),yx=(-0.0080128205128205086393844),yy=(-0.0053418803418803480556609),yz=(-0.0080128205128205069046610),zx=(-0.0801282051282050933327383),zy=(0.0080128205128205069046610),zz=(-0.0854700854700854439904845)}
case (5,1) -> {xx=(0.0534188034188033747384772),xy=(0.0400641025641025327885814),xz=(-0.0160256410256410172787689),yx=(0.0400641025641025327885814),yy=(0.0534188034188033747384772),yz=(-0.0160256410256410172787689),zx=(0.0160256410256410172787689),zy=(0.0160256410256410172787689),zz=(-0.1068376068376068049881056)}
case (5,2) -> {xx=(-0.0053418803418803549945548),xy=(-0.0080128205128205069046610),xz=(-0.0080128205128205069046610),yx=(0.0080128205128205086393844),yy=(-0.0854700854700854439904845),yz=(-0.0801282051282051072105261),zx=(0.0080128205128205069046610),zy=(-0.0801282051282050933327383),zz=(-0.0854700854700854439904845)}
case (5,3) -> {xx=(-0.0587606837606837115184355),xy=(-0.0400641025641025327885814),xz=(-0.0400641025641025397274753),yx=(-0.0400641025641025327885814),yy=(-0.0587606837606837184573294),yz=(-0.0400641025641025397274753),zx=(-0.0400641025641025397274753),zy=(-0.0400641025641025397274753),zz=(-0.0587606837606837184573294)}
case (5,4) -> {xx=(-0.1068376068376067911103178),xy=(0.0160256410256410138093219),xz=(0.0160256410256410172787689),yx=(-0.0160256410256410172787689),yy=(0.0534188034188033955551589),yz=(0.0400641025641025397274753),zx=(-0.0160256410256410172787689),zy=(0.0400641025641025397274753),zz=(0.0534188034188033886162650)}
case (5,5) -> {xx=(0.2350427350427349570960445),xy=(0.0801282051282050933327383),xz=(0.0801282051282050933327383),yx=(0.0801282051282050933327383),yy=(0.2350427350427349570960445),yz=(0.0801282051282050933327383),zx=(0.0801282051282050933327383),zy=(0.0801282051282050933327383),zz=(0.2350427350427349570960445)}
case (5,6) -> {xx=(0.0534188034188033677995833),xy=(-0.0160256410256410138093219),xz=(0.0400641025641025397274753),yx=(0.0160256410256410207482158),yy=(-0.1068376068376068188658934),yz=(0.0160256410256410172787689),zx=(0.0400641025641025397274753),zy=(-0.0160256410256410172787689),zz=(0.0534188034188033816773711)}
case (5,7) -> {xx=(-0.0854700854700854301126967),xy=(-0.0801282051282050933327383),xz=(0.0080128205128205034352140),yx=(-0.0801282051282050933327383),yy=(-0.0854700854700854439904845),yz=(0.0080128205128205086393844),zx=(-0.0080128205128205103741079),zy=(-0.0080128205128205069046610),zz=(-0.0053418803418803419841288)}
case (6,0) -> {xx=(-0.0587606837606837115184355),xy=(0.0400641025641025397274753),xz=(-0.0400641025641025327885814),yx=(0.0400641025641025397274753),yy=(-0.0587606837606837184573294),yz=(0.0400641025641025397274753),zx=(-0.0400641025641025327885814),zy=(0.0400641025641025397274753),zz=(-0.0587606837606837115184355)}
case (6,1) -> {xx=(-0.0053418803418803506577461),xy=(0.0080128205128205086393844),xz=(-0.0080128205128205086393844),yx=(-0.0080128205128205086393844),yy=(-0.0854700854700854301126967),yz=(0.0801282051282050933327383),zx=(0.0080128205128205069046610),zy=(0.0801282051282050933327383),zz=(-0.0854700854700854301126967)}
case (6,2) -> {xx=(0.0534188034188033677995833),xy=(-0.0400641025641025327885814),xz=(-0.0160256410256410207482158),yx=(-0.0400641025641025327885814),yy=(0.0534188034188033677995833),yz=(0.0160256410256410242176628),zx=(0.0160256410256410207482158),zy=(-0.0160256410256410172787689),zz=(-0.1068376068376068466214690)}
case (6,3) -> {xx=(-0.0854700854700854162349088),xy=(-0.0080128205128205086393844),xz=(-0.0801282051282050933327383),yx=(0.0080128205128205086393844),yy=(-0.0053418803418803437188522),yz=(0.0080128205128205121088314),zx=(-0.0801282051282050933327383),zy=(-0.0080128205128205051699375),zz=(-0.0854700854700854439904845)}
case (6,4) -> {xx=(-0.0854700854700854301126967),xy=(0.0801282051282050933327383),xz=(0.0080128205128205121088314),yx=(0.0801282051282050933327383),yy=(-0.0854700854700854439904845),yz=(-0.0080128205128205086393844),zx=(-0.0080128205128205069046610),zy=(0.0080128205128205069046610),zz=(-0.0053418803418803402494053)}
case (6,5) -> {xx=(0.0534188034188033677995833),xy=(0.0160256410256410242176628),xz=(0.0400641025641025466663692),yx=(-0.0160256410256410172787689),yy=(-0.1068376068376068188658934),yz=(-0.0160256410256410207482158),zx=(0.0400641025641025397274753),zy=(0.0160256410256410172787689),zz=(0.0534188034188033955551589)}
case (6,6) -> {xx=(0.2350427350427349293404689),xy=(-0.0801282051282050933327383),xz=(0.0801282051282050933327383),yx=(-0.0801282051282050933327383),yy=(0.2350427350427349570960445),yz=(-0.0801282051282051072105261),zx=(0.0801282051282050933327383),zy=(-0.0801282051282051072105261),zz=(0.2350427350427349848516201)}
case (6,7) -> {xx=(-0.1068376068376067633547422),xy=(-0.0160256410256410172787689),xz=(0.0160256410256410103398750),yx=(0.0160256410256410172787689),yy=(0.0534188034188033955551589),yz=(-0.0400641025641025397274753),zx=(-0.0160256410256410172787689),zy=(-0.0400641025641025397274753),zz=(0.0534188034188033886162650)}
case (7,0) -> {xx=(-0.0053418803418803480556609),xy=(-0.0080128205128205086393844),xz=(0.0080128205128205086393844),yx=(0.0080128205128205086393844),yy=(-0.0854700854700854439904845),yz=(0.0801282051282050933327383),zx=(-0.0080128205128205086393844),zy=(0.0801282051282050933327383),zz=(-0.0854700854700854301126967)}
case (7,1) -> {xx=(-0.0587606837606837115184355),xy=(-0.0400641025641025327885814),xz=(0.0400641025641025397274753),yx=(-0.0400641025641025327885814),yy=(-0.0587606837606837184573294),yz=(0.0400641025641025397274753),zx=(0.0400641025641025397274753),zy=(0.0400641025641025397274753),zz=(-0.0587606837606837184573294)}
case (7,2) -> {xx=(-0.0854700854700854162349088),xy=(0.0080128205128205086393844),xz=(0.0801282051282051072105261),yx=(-0.0080128205128205086393844),yy=(-0.0053418803418803480556609),yz=(0.0080128205128205155782783),zx=(0.0801282051282050794549505),zy=(-0.0080128205128205034352140),zz=(-0.0854700854700854439904845)}
case (7,3) -> {xx=(0.0534188034188033608606894),xy=(0.0400641025641025327885814),xz=(0.0160256410256410207482158),yx=(0.0400641025641025327885814),yy=(0.0534188034188033816773711),yz=(0.0160256410256410242176628),zx=(-0.0160256410256410138093219),zy=(-0.0160256410256410068704280),zz=(-0.1068376068376068188658934)}
case (7,4) -> {xx=(0.0534188034188033677995833),xy=(-0.0160256410256410172787689),xz=(-0.0400641025641025397274753),yx=(0.0160256410256410172787689),yy=(-0.1068376068376068049881056),yz=(-0.0160256410256410172787689),zx=(-0.0400641025641025397274753),zy=(0.0160256410256410172787689),zz=(0.0534188034188033955551589)}
case (7,5) -> {xx=(-0.0854700854700854301126967),xy=(-0.0801282051282050933327383),xz=(-0.0080128205128205121088314),yx=(-0.0801282051282050933327383),yy=(-0.0854700854700854439904845),yz=(-0.0080128205128205086393844),zx=(0.0080128205128205051699375),zy=(0.0080128205128205086393844),zz=(-0.0053418803418803411167670)}
case (7,6) -> {xx=(-0.1068376068376067633547422),xy=(0.0160256410256410172787689),xz=(-0.0160256410256410172787689),yx=(-0.0160256410256410172787689),yy=(0.0534188034188033955551589),yz=(-0.0400641025641025466663692),zx=(0.0160256410256410103398750),zy=(-0.0400641025641025466663692),zz=(0.0534188034188033886162650)}
case (7,7) -> {xx=(0.2350427350427349570960445),xy=(0.0801282051282050933327383),xz=(-0.0801282051282050933327383),yx=(0.0801282051282050933327383),yy=(0.2350427350427350126071957),yz=(-0.0801282051282050933327383),zx=(-0.0801282051282050933327383),zy=(-0.0801282051282050933327383),zz=(0.2350427350427350126071957)}
case _ -> {xx=0,xy=0,xz=0,yx=0,yy=0,yz=0,zx=0,zy=0,zz=0}


let keconst :[24][24]f64 = 
[[0.2350427350427349848516201,-0.0801282051282050933327383,0.0801282051282050933327383,-0.1068376068376068188658934,-0.0160256410256410172787689,0.0160256410256410172787689,-0.0854700854700854301126967,0.0801282051282050933327383,0.0080128205128205121088314,0.0534188034188033816773711,0.0160256410256410172787689,0.0400641025641025327885814,0.0534188034188033816773711,-0.0400641025641025327885814,-0.0160256410256410172787689,-0.0854700854700854301126967,-0.0080128205128205086393844,-0.0801282051282050933327383,-0.0587606837606837115184355,0.0400641025641025397274753,-0.0400641025641025327885814,-0.0053418803418803480556609,0.0080128205128205086393844,-0.0080128205128205051699375],
[-0.0801282051282050933327383,0.2350427350427349848516201,-0.0801282051282050933327383,0.0160256410256410138093219,0.0534188034188033816773711,-0.0400641025641025397274753,0.0801282051282050794549505,-0.0854700854700854301126967,-0.0080128205128205138435549,-0.0160256410256410172787689,-0.1068376068376068049881056,-0.0160256410256410172787689,-0.0400641025641025327885814,0.0534188034188033816773711,0.0160256410256410207482158,0.0080128205128205086393844,-0.0053418803418803480556609,0.0080128205128205103741079,0.0400641025641025397274753,-0.0587606837606837184573294,0.0400641025641025397274753,-0.0080128205128205051699375,-0.0854700854700854439904845,0.0801282051282050933327383],
[0.0801282051282050933327383,-0.0801282051282050933327383,0.2350427350427349570960445,-0.0160256410256410172787689,-0.0400641025641025397274753,0.0534188034188033816773711,-0.0080128205128205034352140,0.0080128205128205034352140,-0.0053418803418803411167670,0.0400641025641025258496875,0.0160256410256410138093219,0.0534188034188033886162650,0.0160256410256410172787689,-0.0160256410256410172787689,-0.1068376068376067911103178,-0.0801282051282050933327383,-0.0080128205128205069046610,-0.0854700854700854439904845,-0.0400641025641025327885814,0.0400641025641025397274753,-0.0587606837606837115184355,0.0080128205128205069046610,0.0801282051282050794549505,-0.0854700854700854301126967],
[-0.1068376068376068049881056,0.0160256410256410138093219,-0.0160256410256410172787689,0.2350427350427349293404689,0.0801282051282050933327383,-0.0801282051282050933327383,0.0534188034188033816773711,-0.0160256410256410172787689,-0.0400641025641025397274753,-0.0854700854700854301126967,-0.0801282051282050794549505,-0.0080128205128205121088314,-0.0854700854700854301126967,0.0080128205128205086393844,0.0801282051282050933327383,0.0534188034188033747384772,0.0400641025641025397274753,0.0160256410256410172787689,-0.0053418803418803506577461,-0.0080128205128205086393844,0.0080128205128205086393844,-0.0587606837606837115184355,-0.0400641025641025327885814,0.0400641025641025397274753],
[-0.0160256410256410172787689,0.0534188034188033816773711,-0.0400641025641025397274753,0.0801282051282050933327383,0.2350427350427349293404689,-0.0801282051282050933327383,0.0160256410256410172787689,-0.1068376068376067911103178,-0.0160256410256410311565567,-0.0801282051282050933327383,-0.0854700854700854301126967,-0.0080128205128205103741079,-0.0080128205128205086393844,-0.0053418803418803454535757,0.0080128205128205086393844,0.0400641025641025397274753,0.0534188034188033747384772,0.0160256410256410172787689,0.0080128205128205086393844,-0.0854700854700854301126967,0.0801282051282050933327383,-0.0400641025641025327885814,-0.0587606837606837184573294,0.0400641025641025397274753],
[0.0160256410256410172787689,-0.0400641025641025327885814,0.0534188034188033816773711,-0.0801282051282050933327383,-0.0801282051282050933327383,0.2350427350427349293404689,-0.0400641025641025327885814,0.0160256410256410103398750,0.0534188034188033816773711,0.0080128205128205017004905,0.0080128205128205034352140,-0.0053418803418803428514905,0.0801282051282050794549505,-0.0080128205128205086393844,-0.0854700854700854301126967,-0.0160256410256410172787689,-0.0160256410256410172787689,-0.1068376068376068049881056,-0.0080128205128205069046610,0.0801282051282050794549505,-0.0854700854700854301126967,0.0400641025641025327885814,0.0400641025641025397274753,-0.0587606837606837184573294],
[-0.0854700854700854301126967,0.0801282051282050794549505,-0.0080128205128205034352140,0.0534188034188033816773711,0.0160256410256410207482158,-0.0400641025641025327885814,0.2350427350427349570960445,-0.0801282051282050933327383,-0.0801282051282051210883139,-0.1068376068376068049881056,-0.0160256410256410172787689,-0.0160256410256410103398750,-0.0587606837606837184573294,0.0400641025641025397274753,0.0400641025641025397274753,-0.0053418803418803558619166,0.0080128205128205069046610,0.0080128205128205086393844,0.0534188034188033677995833,-0.0400641025641025327885814,0.0160256410256410138093219,-0.0854700854700854162349088,-0.0080128205128205051699375,0.0801282051282050794549505],
[0.0801282051282050933327383,-0.0854700854700854301126967,0.0080128205128205051699375,-0.0160256410256410172787689,-0.1068376068376067911103178,0.0160256410256410103398750,-0.0801282051282050933327383,0.2350427350427349570960445,0.0801282051282051210883139,0.0160256410256410172787689,0.0534188034188033816773711,0.0400641025641025397274753,0.0400641025641025397274753,-0.0587606837606837184573294,-0.0400641025641025397274753,-0.0080128205128205086393844,-0.0854700854700854439904845,-0.0801282051282050933327383,-0.0400641025641025327885814,0.0534188034188033816773711,-0.0160256410256410172787689,0.0080128205128205086393844,-0.0053418803418803463209374,-0.0080128205128205017004905],
[0.0080128205128205121088314,-0.0080128205128205103741079,-0.0053418803418803411167670,-0.0400641025641025397274753,-0.0160256410256410242176628,0.0534188034188033816773711,-0.0801282051282051210883139,0.0801282051282051210883139,0.2350427350427350403627713,0.0160256410256410276871097,0.0400641025641025397274753,0.0534188034188033955551589,0.0400641025641025397274753,-0.0400641025641025397274753,-0.0587606837606837253962233,-0.0080128205128205086393844,-0.0801282051282051072105261,-0.0854700854700854439904845,-0.0160256410256410207482158,0.0160256410256410242176628,-0.1068376068376068466214690,0.0801282051282051072105261,0.0080128205128205121088314,-0.0854700854700854439904845],
[0.0534188034188033816773711,-0.0160256410256410172787689,0.0400641025641025327885814,-0.0854700854700854301126967,-0.0801282051282050933327383,0.0080128205128205017004905,-0.1068376068376068049881056,0.0160256410256410172787689,0.0160256410256410276871097,0.2350427350427349848516201,0.0801282051282050933327383,0.0801282051282050933327383,-0.0053418803418803575966400,-0.0080128205128205051699375,-0.0080128205128205069046610,-0.0587606837606837115184355,-0.0400641025641025327885814,-0.0400641025641025397274753,-0.0854700854700854162349088,0.0080128205128205051699375,-0.0801282051282050933327383,0.0534188034188033608606894,0.0400641025641025327885814,-0.0160256410256410103398750],
[0.0160256410256410172787689,-0.1068376068376068188658934,0.0160256410256410172787689,-0.0801282051282050794549505,-0.0854700854700854301126967,0.0080128205128205034352140,-0.0160256410256410172787689,0.0534188034188033816773711,0.0400641025641025397274753,0.0801282051282050933327383,0.2350427350427349848516201,0.0801282051282051072105261,0.0080128205128205086393844,-0.0854700854700854439904845,-0.0801282051282051072105261,-0.0400641025641025327885814,-0.0587606837606837184573294,-0.0400641025641025397274753,-0.0080128205128205051699375,-0.0053418803418803437188522,-0.0080128205128205051699375,0.0400641025641025327885814,0.0534188034188033816773711,-0.0160256410256410103398750],
[0.0400641025641025327885814,-0.0160256410256410172787689,0.0534188034188033886162650,-0.0080128205128205103741079,-0.0080128205128205103741079,-0.0053418803418803437188522,-0.0160256410256410068704280,0.0400641025641025397274753,0.0534188034188034094329467,0.0801282051282050933327383,0.0801282051282051072105261,0.2350427350427349848516201,0.0080128205128205086393844,-0.0801282051282050933327383,-0.0854700854700854439904845,-0.0400641025641025397274753,-0.0400641025641025397274753,-0.0587606837606837184573294,-0.0801282051282050794549505,0.0080128205128205121088314,-0.0854700854700854439904845,0.0160256410256410207482158,0.0160256410256410242176628,-0.1068376068376068188658934],
[0.0534188034188033816773711,-0.0400641025641025327885814,0.0160256410256410172787689,-0.0854700854700854301126967,-0.0080128205128205086393844,0.0801282051282050933327383,-0.0587606837606837184573294,0.0400641025641025397274753,0.0400641025641025397274753,-0.0053418803418803575966400,0.0080128205128205069046610,0.0080128205128205069046610,0.2350427350427349848516201,-0.0801282051282050933327383,-0.0801282051282050933327383,-0.1068376068376067911103178,-0.0160256410256410138093219,-0.0160256410256410207482158,-0.0854700854700854301126967,0.0801282051282050933327383,-0.0080128205128205034352140,0.0534188034188033677995833,0.0160256410256410172787689,-0.0400641025641025327885814],
[-0.0400641025641025327885814,0.0534188034188033816773711,-0.0160256410256410172787689,0.0080128205128205086393844,-0.0053418803418803454535757,-0.0080128205128205086393844,0.0400641025641025397274753,-0.0587606837606837184573294,-0.0400641025641025397274753,-0.0080128205128205086393844,-0.0854700854700854439904845,-0.0801282051282050933327383,-0.0801282051282050933327383,0.2350427350427349848516201,0.0801282051282050933327383,0.0160256410256410172787689,0.0534188034188033955551589,0.0400641025641025466663692,0.0801282051282050933327383,-0.0854700854700854439904845,0.0080128205128205086393844,-0.0160256410256410172787689,-0.1068376068376068188658934,0.0160256410256410172787689],
[-0.0160256410256410172787689,0.0160256410256410207482158,-0.1068376068376067911103178,0.0801282051282050933327383,0.0080128205128205069046610,-0.0854700854700854301126967,0.0400641025641025397274753,-0.0400641025641025397274753,-0.0587606837606837253962233,-0.0080128205128205069046610,-0.0801282051282051072105261,-0.0854700854700854439904845,-0.0801282051282050933327383,0.0801282051282050933327383,0.2350427350427349848516201,0.0160256410256410103398750,0.0400641025641025466663692,0.0534188034188033886162650,0.0080128205128205103741079,-0.0080128205128205086393844,-0.0053418803418803428514905,-0.0400641025641025397274753,-0.0160256410256410172787689,0.0534188034188033955551589],
[-0.0854700854700854301126967,0.0080128205128205069046610,-0.0801282051282050933327383,0.0534188034188033747384772,0.0400641025641025327885814,-0.0160256410256410172787689,-0.0053418803418803549945548,-0.0080128205128205069046610,-0.0080128205128205069046610,-0.0587606837606837115184355,-0.0400641025641025327885814,-0.0400641025641025397274753,-0.1068376068376067911103178,0.0160256410256410138093219,0.0160256410256410172787689,0.2350427350427349570960445,0.0801282051282050933327383,0.0801282051282050933327383,0.0534188034188033677995833,-0.0160256410256410138093219,0.0400641025641025397274753,-0.0854700854700854301126967,-0.0801282051282050933327383,0.0080128205128205034352140],
[-0.0080128205128205086393844,-0.0053418803418803480556609,-0.0080128205128205069046610,0.0400641025641025327885814,0.0534188034188033747384772,-0.0160256410256410172787689,0.0080128205128205086393844,-0.0854700854700854439904845,-0.0801282051282051072105261,-0.0400641025641025327885814,-0.0587606837606837184573294,-0.0400641025641025397274753,-0.0160256410256410172787689,0.0534188034188033955551589,0.0400641025641025397274753,0.0801282051282050933327383,0.2350427350427349570960445,0.0801282051282050933327383,0.0160256410256410207482158,-0.1068376068376068188658934,0.0160256410256410172787689,-0.0801282051282050933327383,-0.0854700854700854439904845,0.0080128205128205086393844],
[-0.0801282051282050933327383,0.0080128205128205069046610,-0.0854700854700854439904845,0.0160256410256410172787689,0.0160256410256410172787689,-0.1068376068376068049881056,0.0080128205128205069046610,-0.0801282051282050933327383,-0.0854700854700854439904845,-0.0400641025641025397274753,-0.0400641025641025397274753,-0.0587606837606837184573294,-0.0160256410256410172787689,0.0400641025641025397274753,0.0534188034188033886162650,0.0801282051282050933327383,0.0801282051282050933327383,0.2350427350427349570960445,0.0400641025641025397274753,-0.0160256410256410172787689,0.0534188034188033816773711,-0.0080128205128205103741079,-0.0080128205128205069046610,-0.0053418803418803419841288],
[-0.0587606837606837115184355,0.0400641025641025397274753,-0.0400641025641025327885814,-0.0053418803418803506577461,0.0080128205128205086393844,-0.0080128205128205086393844,0.0534188034188033677995833,-0.0400641025641025327885814,-0.0160256410256410207482158,-0.0854700854700854162349088,-0.0080128205128205086393844,-0.0801282051282050933327383,-0.0854700854700854301126967,0.0801282051282050933327383,0.0080128205128205121088314,0.0534188034188033677995833,0.0160256410256410242176628,0.0400641025641025466663692,0.2350427350427349293404689,-0.0801282051282050933327383,0.0801282051282050933327383,-0.1068376068376067633547422,-0.0160256410256410172787689,0.0160256410256410103398750],
[0.0400641025641025397274753,-0.0587606837606837184573294,0.0400641025641025397274753,-0.0080128205128205086393844,-0.0854700854700854301126967,0.0801282051282050933327383,-0.0400641025641025327885814,0.0534188034188033677995833,0.0160256410256410242176628,0.0080128205128205086393844,-0.0053418803418803437188522,0.0080128205128205121088314,0.0801282051282050933327383,-0.0854700854700854439904845,-0.0080128205128205086393844,-0.0160256410256410172787689,-0.1068376068376068188658934,-0.0160256410256410207482158,-0.0801282051282050933327383,0.2350427350427349570960445,-0.0801282051282051072105261,0.0160256410256410172787689,0.0534188034188033955551589,-0.0400641025641025397274753],
[-0.0400641025641025327885814,0.0400641025641025397274753,-0.0587606837606837115184355,0.0080128205128205069046610,0.0801282051282050933327383,-0.0854700854700854301126967,0.0160256410256410207482158,-0.0160256410256410172787689,-0.1068376068376068466214690,-0.0801282051282050933327383,-0.0080128205128205051699375,-0.0854700854700854439904845,-0.0080128205128205069046610,0.0080128205128205069046610,-0.0053418803418803402494053,0.0400641025641025397274753,0.0160256410256410172787689,0.0534188034188033955551589,0.0801282051282050933327383,-0.0801282051282051072105261,0.2350427350427349848516201,-0.0160256410256410172787689,-0.0400641025641025397274753,0.0534188034188033886162650],
[-0.0053418803418803480556609,-0.0080128205128205086393844,0.0080128205128205086393844,-0.0587606837606837115184355,-0.0400641025641025327885814,0.0400641025641025397274753,-0.0854700854700854162349088,0.0080128205128205086393844,0.0801282051282051072105261,0.0534188034188033608606894,0.0400641025641025327885814,0.0160256410256410207482158,0.0534188034188033677995833,-0.0160256410256410172787689,-0.0400641025641025397274753,-0.0854700854700854301126967,-0.0801282051282050933327383,-0.0080128205128205121088314,-0.1068376068376067633547422,0.0160256410256410172787689,-0.0160256410256410172787689,0.2350427350427349570960445,0.0801282051282050933327383,-0.0801282051282050933327383],
[0.0080128205128205086393844,-0.0854700854700854439904845,0.0801282051282050933327383,-0.0400641025641025327885814,-0.0587606837606837184573294,0.0400641025641025397274753,-0.0080128205128205086393844,-0.0053418803418803480556609,0.0080128205128205155782783,0.0400641025641025327885814,0.0534188034188033816773711,0.0160256410256410242176628,0.0160256410256410172787689,-0.1068376068376068049881056,-0.0160256410256410172787689,-0.0801282051282050933327383,-0.0854700854700854439904845,-0.0080128205128205086393844,-0.0160256410256410172787689,0.0534188034188033955551589,-0.0400641025641025466663692,0.0801282051282050933327383,0.2350427350427350126071957,-0.0801282051282050933327383],
[-0.0080128205128205086393844,0.0801282051282050933327383,-0.0854700854700854301126967,0.0400641025641025397274753,0.0400641025641025397274753,-0.0587606837606837184573294,0.0801282051282050794549505,-0.0080128205128205034352140,-0.0854700854700854439904845,-0.0160256410256410138093219,-0.0160256410256410068704280,-0.1068376068376068188658934,-0.0400641025641025397274753,0.0160256410256410172787689,0.0534188034188033955551589,0.0080128205128205051699375,0.0080128205128205086393844,-0.0053418803418803411167670,0.0160256410256410103398750,-0.0400641025641025466663692,0.0534188034188033886162650,-0.0801282051282050933327383,-0.0801282051282050933327383,0.2350427350427350126071957]
]
