{
  description = "Replicating and Comparing the PIP-ECHO Paper";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.05";
    openms.url = "github:pjones/nix-openms";
  };

  outputs = { self, nixpkgs, ... }:
    let
      supportedSystems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];

      # Function to generate a set based on supported systems:
      forAllSystems = f:
        nixpkgs.lib.genAttrs supportedSystems (system: f system);

      # Attribute set of nixpkgs for each system:
      nixpkgsFor = forAllSystems (system:
        import nixpkgs { inherit system; });

      # Packages needed to build the wiki:
      buildInputs = pkgs: with pkgs; [
        (python3.withPackages (ps: with ps; [
          pygments
        ]))

        curl
        texliveFull

        # OpenMS
        self.inputs.openms.packages.${pkgs.system}.openms
      ];
    in
    {
      devShells = forAllSystems (system:
        let pkgs = nixpkgsFor.${system}; in
        {
          default = pkgs.mkShell {
            NIX_PATH = "nixpkgs=${pkgs.path}";
            buildInputs = buildInputs pkgs;
          };
        });
    };
}
