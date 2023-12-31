module Numerics where
    import qualified Data.Map as M
    import qualified Data.Set as S

    data Element = 
        Carbon | 
        Hydrogen | 
        Oxygen |
        Nitrogen | 
        Helium |
        Fluorine |
        Sulfur |
        Chlorine |
        Barium | 
        Copper |
        Yttrium |
        Mercury  | 
        Phosphorus |
        Iron |
        Aluminium |
        Sodium
        deriving (Eq, Ord, Show)

    data Compound = 
        Hydrogen2 |
        CarbonDioxide | 
        Helium2 |
        SulfurHexafluoride |
        Fluorine2 |
        ChloralHydrate |
        Fluoxymesterone | 
        Formaldehyde | 
        Glucose | 
        AceticAcid | 
        Water |
        Ammonia |
        Oxygen2 | 
        Iron3Sulfide |
        AmmoniumPerchlorate |
        AluminiumOxide |
        AluminiumChloride |
        NitricOxide |
        SodiumBicarbonate |
        CitricAcid | 
        BariumPeroxide |
        Bornite | 
        HydrogenCyanide |
        Methane
        deriving (Eq, Show)

    type EmpiricalFormula = [(Element, Integer)]
    type EmpiricalFormulaMultiple = (EmpiricalFormula, Integer)

    data Reaction = Reaction [EmpiricalFormulaMultiple] [EmpiricalFormulaMultiple]
    
    relativeMass :: Element -> Double
    relativeMass Carbon = 12.011
    relativeMass Hydrogen = 1.0078
    relativeMass Oxygen = 15.999
    relativeMass Nitrogen = 14.007
    relativeMass Helium = 4.0026
    relativeMass Fluorine = 18.998
    relativeMass Sulfur = 32.065
    relativeMass Chlorine = 35.453
    relativeMass Barium = 137.33 
    relativeMass Copper = 63.546
    relativeMass Yttrium = 88.906
    relativeMass Mercury = 200.59
    relativeMass Phosphorus = 30.974
    relativeMass Iron = 55.845
    relativeMass Aluminium = 26.982
    relativeMass Sodium = 22.990 

    empiricalFormula :: Compound -> EmpiricalFormula
    empiricalFormula Hydrogen2 = [(Hydrogen, 2)]
    empiricalFormula CarbonDioxide = [(Carbon, 1), (Oxygen, 2)]
    empiricalFormula Helium2 = [(Helium, 2)]
    empiricalFormula SulfurHexafluoride = [(Sulfur, 1), (Fluorine, 6)]
    empiricalFormula Fluorine2 = [(Fluorine, 2)]
    empiricalFormula ChloralHydrate = [(Carbon, 2), (Hydrogen, 3), (Chlorine, 3), (Oxygen, 2)]
    empiricalFormula Fluoxymesterone = [(Carbon, 20), (Hydrogen, 29), (Fluorine, 1), (Oxygen, 3)]
    empiricalFormula Formaldehyde = [(Carbon, 1), (Hydrogen, 1), (Oxygen, 1)]
    empiricalFormula Glucose = [(Carbon, 6), (Hydrogen, 12), (Oxygen, 6)]
    empiricalFormula AceticAcid = [(Hydrogen, 4), (Carbon, 2), (Oxygen, 6)]
    empiricalFormula Water = [(Hydrogen, 2), (Oxygen, 1)]
    empiricalFormula Ammonia = [(Nitrogen, 1), (Hydrogen, 3)]
    empiricalFormula Oxygen2 = [(Oxygen, 2)]
    empiricalFormula Iron3Sulfide = [(Iron, 2), (Sulfur, 3)]
    empiricalFormula AmmoniumPerchlorate = [(Nitrogen, 1), (Hydrogen, 4), (Chlorine, 1), (Oxygen, 4)]
    empiricalFormula AluminiumOxide = [(Aluminium, 2), (Oxygen, 3)]
    empiricalFormula AluminiumChloride = [(Aluminium, 1), (Chlorine, 3)]
    empiricalFormula NitricOxide = [(Nitrogen, 1), (Oxygen, 1)]
    empiricalFormula SodiumBicarbonate = [(Sodium, 1), (Hydrogen, 1), (Carbon, 1), (Oxygen, 3)]
    empiricalFormula CitricAcid = [(Carbon, 6), (Hydrogen, 8), (Oxygen, 7)]
    empiricalFormula BariumPeroxide = [(Barium, 1), (Oxygen, 2)]
    empiricalFormula Bornite = [(Copper, 3), (Iron, 1), (Sulfur, 3)]
    empiricalFormula HydrogenCyanide = [(Hydrogen, 1), (Carbon, 1), (Nitrogen, 1)]
    empiricalFormula Methane = [(Carbon, 1), (Hydrogen, 4)]

    empForMulToEmpFor :: EmpiricalFormulaMultiple -> EmpiricalFormula
    empForMulToEmpFor (emf, n) = map (\(e, m) -> (e, m*n)) emf

    mapToEmpiricalFormula :: M.Map Element Integer -> EmpiricalFormula
    mapToEmpiricalFormula = M.foldrWithKey (\x y xs -> (x,y):xs) [] 

    empiricalFormulaToMap :: EmpiricalFormula -> M.Map Element Integer
    empiricalFormulaToMap = foldr (\(e, n) xs -> M.insert e n xs) M.empty 

    addEmpiricalFormula :: EmpiricalFormula -> EmpiricalFormula -> EmpiricalFormula
    addEmpiricalFormula emf1 emf2 = mapToEmpiricalFormula $ foldr (\(e, n) xs -> M.insert e (M.findWithDefault 0 e xs + n) xs) M.empty (emf1 ++ emf2)

    glucoseCombustionReaction = Reaction [(empiricalFormula Glucose, 1), (empiricalFormula Oxygen2, 6)] [(empiricalFormula CarbonDioxide, 6), (empiricalFormula Water, 6)]


    checkReaction :: Reaction -> Bool
    checkReaction (Reaction left right) = (S.isSubsetOf left' right') && (S.isSubsetOf right' left')
                                            where 
                                                f emfm = S.fromList . M.toList . empiricalFormulaToMap $ foldr (\x xs -> addEmpiricalFormula x xs) [] (map empForMulToEmpFor emfm)
                                                left' = f left
                                                right' = f right

    avogadroNumber :: Double
    avogadroNumber = 6.02214076e23

    data MassUnit = 
        MiliGram Double | 
        Gram Double |
        KiloGram Double |
        MetricTon Double
        deriving (Eq, Show)

    data AmUnit =
        Mol Double |
        TheNumber Double
        deriving (Eq, Show)

    data VolumeUnit = 
        MiliLitter Double |
        Litter Double 
        deriving (Eq, Ord, Show)

    sumMass :: EmpiricalFormula -> Double
    sumMass = sum . map (\(e,n) -> relativeMass e * (fromInteger n)) 

    massToGram :: MassUnit -> Double
    massToGram (MiliGram x) = x / 1000
    massToGram (Gram x) = x
    massToGram (KiloGram x) = 1000 * x
    massToGram (MetricTon x) = 1000 * 1000 * x

    massToMol :: MassUnit -> Double -> Double
    massToMol m = (/) (massToGram m)

    amountToMol :: AmUnit -> Double
    amountToMol (Mol x) = x 
    amountToMol (TheNumber x) = x / avogadroNumber

    amountToMass :: AmUnit -> Double -> Double
    amountToMass am m = m * amountToMol am

    elementToMol :: MassUnit -> Element -> Double
    elementToMol m e = massToGram m / relativeMass e

    elementsToMols :: [MassUnit] -> [Element] -> [Double]
    elementsToMols ms es = zipWith elementToMol ms es

    calcEmpricalFormula :: [AmUnit] -> [Element] -> EmpiricalFormula
    calcEmpricalFormula ams es = zip es (map (\x -> round (x / (minimum xs))) xs)
                                    where
                                        xs = map amountToMol ams
    
    calcEmpricalFormulaFromMass :: [MassUnit] -> [Element] -> EmpiricalFormula
    calcEmpricalFormulaFromMass ms es = calcEmpricalFormula (fmap Mol (elementsToMols ms es)) es

    molForToNumAtoms :: EmpiricalFormula -> Integer
    molForToNumAtoms = sum . map snd

    empiricalMol :: EmpiricalFormula -> MassUnit -> Double
    empiricalMol e m = massToMol m (sumMass e)

    compoundMolarMass :: Compound -> Double
    compoundMolarMass = sumMass . empiricalFormula

    compoundMass :: AmUnit -> Compound -> Double
    compoundMass am c = amountToMol am * compoundMolarMass c

    compoundMol :: MassUnit -> Compound -> Double
    compoundMol m c = empiricalMol (empiricalFormula c) m

    compoundMols :: MassUnit -> Compound -> [(Element, Double)]
    compoundMols m c = map (\(e, n) -> (e, (fromInteger n) * x)) (empiricalFormula c)
                        where 
                            x = compoundMol m c

    percentComposition :: EmpiricalFormula -> [(Element, Double)]
    percentComposition mf = map (\(e, n) -> (e, relativeMass e * (fromInteger n) / total)) mf
                                where 
                                    total = sumMass mf

    elementsMass :: EmpiricalFormula -> MassUnit -> [(Element, Double)]
    elementsMass ef m = map (\(e, x) -> (e, x * (massToGram m))) (percentComposition ef)