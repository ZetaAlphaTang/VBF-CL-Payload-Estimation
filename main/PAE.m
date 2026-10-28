function a = PAE(specificForce,thrust,bodyZ,mU,mL,g)
fCable=mU*specificForce-bodyZ.*thrust;
a=g-fCable/mL;
end
