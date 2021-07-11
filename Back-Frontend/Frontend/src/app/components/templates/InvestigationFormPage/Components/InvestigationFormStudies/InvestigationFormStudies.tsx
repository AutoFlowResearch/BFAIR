import React, { useContext } from 'react';
import InvestigationFormStudy from '../InvestigationFormStudy';
import addIcon from '../../../../../../assets/images/add-icon.svg';
import { addStudy, removeStudy } from '../../../../../shared/redux/actions';
import { connect } from 'react-redux';

const InvestigationFormStudies = (props) => {
  const { investigation, addStudy, removeStudy } = props;
  return (
    <React.Fragment>
      {investigation.studies.map((study, index) => (
        <InvestigationFormStudy
          study={study}
          key={index}
          deleteStudy={() => removeStudy(index)}
        />
      ))}
      <div className="section__add">
        <button className="section__add__button" onClick={() => addStudy()}>
          <img alt="" src={addIcon} className="section__add__button__icon" />
          Add New Study
        </button>
      </div>
    </React.Fragment>
  );
};

const mapStateToProps = (state) => ({
  investigation: state,
});

const mapDispatchToProps = (dispatch) => ({
  removeStudy: (index) => dispatch(removeStudy(index)),
  addStudy: () => dispatch(addStudy()),
});

export default connect(
  mapStateToProps,
  mapDispatchToProps,
)(InvestigationFormStudies);
