{
  description = "Replicating and Comparing the PIP-ECHO Paper";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-26.05";
    proteomics.url = "github:pjones/proteomics.nix/openms-next";
  };

  outputs =
    { self, nixpkgs, ... }:
    let
      supportedSystems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];

      # Function to generate a set based on supported systems:
      forAllSystems = f: nixpkgs.lib.genAttrs supportedSystems (system: f system);

      # Attribute set of nixpkgs for each system:
      nixpkgsFor = forAllSystems (system: import nixpkgs { inherit system; });

      # Packages needed to build the wiki:
      buildInputs =
        pkgs: with pkgs; [
          (python3.withPackages (
            ps: with ps; [
              pygments
            ]
          ))

          curl
          jre
          texliveFull

          self.inputs.proteomics.packages.${pkgs.system}.comet
          self.inputs.proteomics.packages.${pkgs.system}.flashlfq
          self.inputs.proteomics.packages.${pkgs.system}.metamorpheus
          self.inputs.proteomics.packages.${pkgs.system}.msgfplus
          self.inputs.proteomics.packages.${pkgs.system}.openms
          self.inputs.proteomics.packages.${pkgs.system}.percolator
          self.inputs.proteomics.packages.${pkgs.system}.thermorawfp

        ];
    in
    {
      devShells = forAllSystems (
        system:
        let
          pkgs = nixpkgsFor.${system};
          lib = pkgs.lib;
          gccVer = lib.concatStringsSep "." (lib.take 3 (lib.splitVersion pkgs.libgcc.version));
        in
        {
          default = pkgs.mkShell {
            NIX_PATH = "nixpkgs=${pkgs.path}";
            buildInputs = buildInputs pkgs;

            # Help GDB find the STL pretty printer:
            PYTHON_LIBSTDCXX = "${pkgs.libgcc.lib}/share/gcc-${gccVer}/python";
          };
        }
      );
    };
}
