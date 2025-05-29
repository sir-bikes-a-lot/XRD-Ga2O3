function [ params ] = binaries( binary )
% BINARIES
% Get the parameters for a specified binary compound. Return a structure of
% the following form:
%
% params = struct(  'a', a,...          a lattice parameter in Angstroms
%                   'b', b,...          b lattice parameter in Angstroms
%                   'c', c,...          c lattice parameter in Angstroms
%                   'beta', beta,...    beta lattice parameter in degrees
%                   'C11', C11,...      Elastic constant in GPa
%                   'C12', C12,...      Elastic constant in GPa
%                   'C13', C13,...      Elastic constant in GPa
%                   'C15', C15,...      Elastic constant in GPa
%                   'C22', C22,...      Elastic constant in GPa
%                   'C23', C23,...      Elastic constant in GPa
%                   'C25', C25,...      Elastic constant in GPa
%                   'C33', C33,...      Elastic constant in GPa
%                   'C35', C35,...      Elastic constant in GPa
%                   'C55', C55);        Elastic constant in GPa
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(binary)
    case('Ga2O3')
%         params = struct('a', 12.21,...      % Ref needed here
%                         'b', 3.04,...       % Ref needed here
%                         'c', 5.81,...       % Ref needed here
%                         'beta', 103.87,...  % Ref needed here
%         params = struct('a', 12.23,...      % E. Hettiaratchy et al, "Quantitative x-ray diffraction analysis of strain and interdiffusion in β-Ga2O3 superlattices of μ-Fe2O3 and β-(AlxGa1−x)2O3", J. Vac. Sci. Technol. A 40, 062708 (2022).
%                         'b', 3.04,...       % E. Hettiaratchy et al, "Quantitative x-ray diffraction analysis of strain and interdiffusion in β-Ga2O3 superlattices of μ-Fe2O3 and β-(AlxGa1−x)2O3", J. Vac. Sci. Technol. A 40, 062708 (2022).
%                         'c', 5.80,...       % E. Hettiaratchy et al, "Quantitative x-ray diffraction analysis of strain and interdiffusion in β-Ga2O3 superlattices of μ-Fe2O3 and β-(AlxGa1−x)2O3", J. Vac. Sci. Technol. A 40, 062708 (2022).
%                         'beta', 103.87,...  % E. Hettiaratchy et al, "Quantitative x-ray diffraction analysis of strain and interdiffusion in β-Ga2O3 superlattices of μ-Fe2O3 and β-(AlxGa1−x)2O3", J. Vac. Sci. Technol. A 40, 062708 (2022).
%        BEST PARAMETERS HERE:
       params = struct('a', 12.214,...     % J. Ahman, G. Svensson, and J. Albertsson, “A Reinvestigation of Beta-Gallium Oxide", Acta Cryst. C52, 1336-1338 (1996).
                        'b', 3.0371,...     % J. Ahman, G. Svensson, and J. Albertsson, “A Reinvestigation of Beta-Gallium Oxide", Acta Cryst. C52, 1336-1338 (1996).
                        'c', 5.7981,...     % J. Ahman, G. Svensson, and J. Albertsson, “A Reinvestigation of Beta-Gallium Oxide", Acta Cryst. C52, 1336-1338 (1996).
                        'beta', 103.83,...  % J. Ahman, G. Svensson, and J. Albertsson, “A Reinvestigation of Beta-Gallium Oxide", Acta Cryst. C52, 1336-1338 (1996).
                        'C11', 214.3,...      % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C12', 110.3,...      % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C13', 120.0,...      % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C15', -19.7,...      % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C22', 330.0,...      % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C23', 66.9,...       % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C25', 11.9,...       % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C33', 324.8,...      % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C35', 7.5,...        % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C55', 68.9);         % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
%                         'C11', 237,...      % A. Mauze, T. Itoh, Y. Zhang, E. Deagueros, F. Wu, J. Speck, J. Appl. Phys. 132, 115302 (2022)
%                         'C12', 125,...      % A. Mauze, T. Itoh, Y. Zhang, E. Deagueros, F. Wu, J. Speck, J. Appl. Phys. 132, 115302 (2022)
%                         'C13', 147,...      % A. Mauze, T. Itoh, Y. Zhang, E. Deagueros, F. Wu, J. Speck, J. Appl. Phys. 132, 115302 (2022)
%                         'C15', -18,...      % A. Mauze, T. Itoh, Y. Zhang, E. Deagueros, F. Wu, J. Speck, J. Appl. Phys. 132, 115302 (2022)
%                         'C22', 354,...      % A. Mauze, T. Itoh, Y. Zhang, E. Deagueros, F. Wu, J. Speck, J. Appl. Phys. 132, 115302 (2022)
%                         'C23', 95,...       % A. Mauze, T. Itoh, Y. Zhang, E. Deagueros, F. Wu, J. Speck, J. Appl. Phys. 132, 115302 (2022)
%                         'C25', 11,...       % A. Mauze, T. Itoh, Y. Zhang, E. Deagueros, F. Wu, J. Speck, J. Appl. Phys. 132, 115302 (2022)
%                         'C33', 257,...      % A. Mauze, T. Itoh, Y. Zhang, E. Deagueros, F. Wu, J. Speck, J. Appl. Phys. 132, 115302 (2022)
%                         'C35', 6,...        % A. Mauze, T. Itoh, Y. Zhang, E. Deagueros, F. Wu, J. Speck, J. Appl. Phys. 132, 115302 (2022)
%                         'C55', 67);         % A. Mauze, T. Itoh, Y. Zhang, E. Deagueros, F. Wu, J. Speck, J. Appl. Phys. 132, 115302 (2022)
    case('Al2O3')
%         params = struct('a', 11.79,...      % C. Kranert, M. Jenderka, J. Lenzner, M. Lorenz, H. von Wenckstern, R. Schmidt-Grund, and M. Grundmann, J. Appl. Phys. 117, 125703 (2015)
%                         'b', 2.91,...       % C. Kranert, M. Jenderka, J. Lenzner, M. Lorenz, H. von Wenckstern, R. Schmidt-Grund, and M. Grundmann, J. Appl. Phys. 117, 125703 (2015)
%                         'c', 5.64,...       % C. Kranert, M. Jenderka, J. Lenzner, M. Lorenz, H. von Wenckstern, R. Schmidt-Grund, and M. Grundmann, J. Appl. Phys. 117, 125703 (2015)
%                         'beta', 104.14,...  % C. Kranert, M. Jenderka, J. Lenzner, M. Lorenz, H. von Wenckstern, R. Schmidt-Grund, and M. Grundmann, J. Appl. Phys. 117, 125703 (2015)
%         params = struct('a', 12.23,...      % E. Hettiaratchy et al, "Quantitative x-ray diffraction analysis of strain and interdiffusion in β-Ga2O3 superlattices of μ-Fe2O3 and β-(AlxGa1−x)2O3", J. Vac. Sci. Technol. A 40, 062708 (2022).
%                         'b', 2.849,...      % E. Hettiaratchy et al, "Quantitative x-ray diffraction analysis of strain and interdiffusion in β-Ga2O3 superlattices of μ-Fe2O3 and β-(AlxGa1−x)2O3", J. Vac. Sci. Technol. A 40, 062708 (2022).
%                         'c', 5.80,...       % E. Hettiaratchy et al, "Quantitative x-ray diffraction analysis of strain and interdiffusion in β-Ga2O3 superlattices of μ-Fe2O3 and β-(AlxGa1−x)2O3", J. Vac. Sci. Technol. A 40, 062708 (2022).
%                         'beta', 103.87,...  % E. Hettiaratchy et al, "Quantitative x-ray diffraction analysis of strain and interdiffusion in β-Ga2O3 superlattices of μ-Fe2O3 and β-(AlxGa1−x)2O3", J. Vac. Sci. Technol. A 40, 062708 (2022).
%         BEST PARAMETERS HERE
        params = struct('a', 11.855,...     % R.-S. Zhou and R. L. Snyder, “Structures and Transformation Mechanisms of the eta, gamma, and theta Transition Aluminas”, Acta Cryst. B47, 617-630 (1991).
                        'b', 2.904,...      % R.-S. Zhou and R. L. Snyder, “Structures and Transformation Mechanisms of the eta, gamma, and theta Transition Aluminas”, Acta Cryst. B47, 617-630 (1991).
                        'c', 5.622,...      % R.-S. Zhou and R. L. Snyder, “Structures and Transformation Mechanisms of the eta, gamma, and theta Transition Aluminas”, Acta Cryst. B47, 617-630 (1991).
                        'beta', 103.84,...  % R.-S. Zhou and R. L. Snyder, “Structures and Transformation Mechanisms of the eta, gamma, and theta Transition Aluminas”, Acta Cryst. B47, 617-630 (1991).
                        'C11', 264.8,...      % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C12', 111.3,...      % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C13', 129.0,...      % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C15', -27.3,...      % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C22', 384.6,...      % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C23', 58.9,...       % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C25', 14.3,...       % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C33', 399.6,...      % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C35', 14.6,...        % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
                        'C55', 101.8);         % R. Korlacki et al, "Strain and Composition Dependencies of the Near-Band-Gap Optical Transitions in Monoclinic (AlxGa1−x)2O3 Alloys with Coherent Biaxial In-Plane Strain on Ga2O3(010)", Phys. Rev. Appl. 18, 064019 (2022).
%                         'C11', 237,...      % Unknown, assume same as Ga2O3
%                         'C12', 125,...      % Unknown, assume same as Ga2O3
%                         'C13', 147,...      % Unknown, assume same as Ga2O3
%                         'C15', -18,...      % Unknown, assume same as Ga2O3
%                         'C22', 354,...      % Unknown, assume same as Ga2O3
%                         'C23', 95,...       % Unknown, assume same as Ga2O3
%                         'C25', 11,...       % Unknown, assume same as Ga2O3
%                         'C33', 257,...      % Unknown, assume same as Ga2O3
%                         'C35', 6,...        % Unknown, assume same as Ga2O3
%                         'C55', 67);         % Unknown, assume same as Ga2O3
    case('In2O3')
        params = struct('a', 13.084,...     % X. Liu and C.-K. Tan, “Electronic properties of monoclinic (InxGa1-x)2O3 alloys by first-principle, AIP Advances 9, 035318 (2019).
                        'b', 3.269,...      % X. Liu and C.-K. Tan, “Electronic properties of monoclinic (InxGa1-x)2O3 alloys by first-principle, AIP Advances 9, 035318 (2019).
                        'c', 6.397,...      % X. Liu and C.-K. Tan, “Electronic properties of monoclinic (InxGa1-x)2O3 alloys by first-principle, AIP Advances 9, 035318 (2019).
                        'beta', 103.83,...  % X. Liu and C.-K. Tan, “Electronic properties of monoclinic (InxGa1-x)2O3 alloys by first-principle, AIP Advances 9, 035318 (2019).
                        'C11', 237,...      % Unknown, assume same as Ga2O3
                        'C12', 125,...      % Unknown, assume same as Ga2O3
                        'C13', 147,...      % Unknown, assume same as Ga2O3
                        'C15', -18,...      % Unknown, assume same as Ga2O3
                        'C22', 354,...      % Unknown, assume same as Ga2O3
                        'C23', 95,...       % Unknown, assume same as Ga2O3
                        'C25', 11,...       % Unknown, assume same as Ga2O3
                        'C33', 54,...       % Unknown, assume same as Ga2O3
                        'C35', 0,...        % Unknown, assume same as Ga2O3
                        'C55', 67);         % Unknown, assume same as Ga2O3
    otherwise
%         error('ERRCODE010: Binary compound not found!');

end

