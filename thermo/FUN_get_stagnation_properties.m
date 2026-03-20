function [p0, T0] = FUN_get_stagnation_properties(T, p, V)
% Frozen-composition stagnation properties for thermally perfect AIR.
% Inputs: T [K], p [Pa], V [m/s] (scalars)
% Outputs: p0 [Pa], T0 [K]

  % constants
  Ru   = 8.31446261815324; 
  Mair = 28.9651159e-3;         
  R    = Ru / Mair;

  h_of_T  = @(Tq) FUN_get_prop_NASA9(Tq,'h');
  s0_of_T = @(Tq) FUN_get_prop_NASA9(Tq,'s');
  cp_of_T = @(Tq) FUN_get_prop_NASA9(Tq,'cp');   % derivative dh/dT

  % total enthalpy and static entropy
  h0 = h_of_T(T) + 0.5*V^2;
  s_static = s0_of_T(T) - R*log(p);  % using same 1-Pa reference on both sides

  % Fast, robust inversion: h(T0) = h0 using safeguarded Newton
  if V==0
      T0 = T;
  else
      % seed from energy: ΔT ≈ (V^2/2)/cp(T)
      Tseed = T + 0.5*V^2/max(cp_of_T(T), 1.0);
      % bracket [Tlo, Thi] so that h(Tlo) <= h0 <= h(Thi)
      Tlo = max(50.0, min(T, Tseed));
      hlo = h_of_T(Tlo);
      if hlo > h0, Tlo = T; hlo = h_of_T(Tlo); end

      Thi = max(Tseed, T)+10.0;
      hhi = h_of_T(Thi);
      grow = 0;
      while hhi < h0 && grow < 100
          Thi = min(2.0*Thi, 3.0e4);   % cap to avoid silly temps
          hhi = h_of_T(Thi);
          grow = grow+1;
      end
      % Safeguarded Newton within [Tlo,Thi]
      Tn = min(max(Tseed, Tlo), Thi);
      for k=1:20
          hn = h_of_T(Tn);
          fn = hn - h0;
          if abs(fn) <= 1e-6*max(1.0, h0), break; end

          cpn = cp_of_T(Tn);                  % dh/dT
          % Newton step
          if cpn > 1e-12
              Tnew = Tn - fn/cpn;
          else
              Tnew = (Tlo+Thi)/2;             % degenerate derivative → bisect
          end

          % keep inside bracket; if outside, bisect
          if Tnew <= Tlo || Tnew >= Thi
              Tnew = 0.5*(Tlo+Thi);
          end

          % shrink bracket based on sign
          if fn > 0
              Thi = Tn;
          else
              Tlo = Tn;
          end
          Tn = Tnew;

          % small bracket → done
          if (Thi - Tlo) <= 1e-6*max(1.0, Tn), break; end
      end
      T0 = Tn;
  end

  % stagnation pressure (frozen isentropic from (T,p) to (T0,p0))
  p0 = exp( (s0_of_T(T0) - s_static)/R );
end
