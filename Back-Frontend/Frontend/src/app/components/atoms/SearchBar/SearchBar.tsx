import React, { KeyboardEvent, useState } from 'react';
import './SearchBar.style.scss';
import searchIcon from '../../../../assets/images/search-icon.svg';

const SearchBar = (props) => {
  const { onSearch = () => {} } = props;
  const [searchTerm, setSearchTerm] = useState('');

  const handleSearchInput = (event) => {
    const searchTerm = event.target.value;
    setSearchTerm(searchTerm);
    onSearch(searchTerm);
  };

  return (
    <div className='search-bar'>
      <img alt='' src={searchIcon} className='search-icon' />
      <input
        placeholder='Search...'
        value={searchTerm}
        onChange={handleSearchInput}
      />
    </div>
  );
};

export default SearchBar;
