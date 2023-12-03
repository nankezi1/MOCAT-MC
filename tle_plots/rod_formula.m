function R=rod_formula(e,w)

R = eye(3,3) + sin(w)*cpm(e) + (1-cos(w))*cpm(e)^2;